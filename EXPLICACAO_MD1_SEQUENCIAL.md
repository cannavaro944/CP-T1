# Explicação Técnica — Etapa 1: Otimizações Sequenciais em `src/MD1.cpp`

> **Âmbito desta etapa:** Otimizações estritamente **sequenciais** (sem `#include <omp.h>` e sem diretivas `#pragma omp`), aplicadas sobre o código base original do professor em `src/MD1.cpp`. A interface de entrada/saída (`inputdata.txt`, `*_traj.xyz`, `*_output.txt`, `*_average.txt`) e a física da simulação (potencial de Lennard-Jones, integrador *Velocity Verlet*, paredes elásticas, $N = 5000$, `NumTime = 200` para Árgon) mantêm-se inalteradas.

---

## 1. Resumo das Alterações por Função

| Função | Linhas em `src/MD1.cpp` | Complexidade | Otimizações aplicadas |
|---|---:|---|---|
| `Potential()` | 456–483 | $O(N^2/2)$ | Eliminação de `sqrt()` e `pow()`; cálculo algébrico via `inv_r2` e `inv_r6`; ciclo triangular (`j = i + 1 .. N - 1`); constante `8.0 * epsilon` fora do ciclo; `rxi, ryi, rzi` e `local_pot` em registos. |
| `computeAccelerations()` | 489–545 | $O(N^2/2)$ | Substituição de `pow(rSqd, -7)` e `pow(rSqd, -4)` por `inv_r2`, `inv_r4`, `inv_r8`, `inv_r14`; constante `c24 = 24.0` fora do ciclo; acumuladores locais `axi, ayi, azi` e coordenadas `rxi, ryi, rzi` em registos; reutilização de `fx, fy, fz` e remoção dos ciclos `k = 0..2`. |
| `VelocityVerlet()` | 548–593 | $O(N)$ (+ forças) | Pré-cálculo de `half_dt`, `half_dt2`, `wall_factor` e `inv_area = 1.0 / (6.0 * L * L)`; leitura única de `a[i][j]` no 1.º ciclo; fusão do 2.º meio-passo de velocidade com verificação de paredes elásticas e acumulação de `psum`. |
| `Kinetic()` | 440–452 | $O(N)$ | Pré-cálculo de `half_m = 0.5 * m` fora do ciclo; remoção do ciclo `j = 0..2` e de $N$ multiplicações/divisões repetidas por passo. |
| `MeanSquaredVelocity()` | 421–437 | $O(N)$ | Pré-cálculo de `inv_N = 1.0 / N` e carregamento único de `v[i][0..2]` em escalares locais. |
| `initializeVelocities()` | 596–648 | $O(N)$ | Pré-cálculo de `inv_Nm = 1.0 / (N * m)`; fusão da subtração de `vCM` com o cálculo de `vSqdSum` numa única passagem sobre `v`. |
| `main()` | 282–348 | $O(\text{NumTime})$ | Pré-cálculo das constantes invariantes `VolSI`, `dt_SI`, `temp_scale`, `gc_scale` e `Z_scale` fora do ciclo temporal. |

---

## 2. Detalhe Técnico e Justificação Física/Computacional

### 2.1. `Potential()` — Linhas 456–483 *(Hotspot principal: ~65,4% na versão original)*

**O que mudou:**
- **Antes:** percorria todos os pares $(i, j)$ com $i, j \in [0, N-1]$ e condição `if (j != i)`, calculava `r2` num ciclo `k = 0..2`, chamava `sqrt(r2)`, dividia `sigma / rnorm`, invocava `pow(quot, 12.)` e `pow(quot, 6.)`, e multiplicava por `4 * epsilon` dentro do ciclo interno.
- **Depois:**
  - Percorre apenas pares únicos `i = 0 .. N - 2` e `j = i + 1 .. N - 1`, removendo o ramo condicional `if (j != i)`.
  - Carrega `rxi, ryi, rzi` fora do ciclo em `j` e acumula o somatório do átomo `i` em `local_pot`.
  - Substitui `sqrt` e `pow` por:
    ```cpp
    double r2 = dx * dx + dy * dy + dz * dz;
    double inv_r2 = 1.0 / r2;
    double inv_r6 = inv_r2 * inv_r2 * inv_r2;
    local_pot += (inv_r6 * inv_r6 - inv_r6);
    ```
  - Pré-calcula `const double pot_factor = 8.0 * epsilon;` fora dos ciclos e multiplica o total `Pot` uma única vez no `return`.

**Motivo físico e computacional:**
1. **Simetria de pares (física):** o potencial de Lennard-Jones depende apenas da distância euclidiana $r_{ij} = \|\mathbf{r}_i - \mathbf{r}_j\| = r_{ji}$, logo $V(r_{ij}) = V(r_{ji})$. Como o código original somava ambos os pares $(i, j)$ e $(j, i)$, somar apenas para $i < j$ e multiplicar no fim por $2.0 \times 4.0\,\varepsilon = 8.0\,\varepsilon$ preserva exatamente o mesmo valor físico cortando o número de pares para metade: de $N(N-1) = 24\,995\,000$ para $N(N-1)/2 = 12\,497\,500$ por passo temporal.
2. **Eliminação de `sqrt` e `pow` (computacional):** em unidades naturais ($\sigma = 1$), $(\sigma / \sqrt{r^2})^6 = r^{-6} = (r^{-2})^3$ e $(\sigma / \sqrt{r^2})^{12} = (r^{-6})^2$. Chamadas genéricas à `libm` (`pow`, `sqrt`) têm latência elevada (dezenas de ciclos) e impedem auto-vetorização eficiente; substituí-las por 1 divisão (`1.0 / r2`) e 3 multiplicações em vírgula flutuante reduz drasticamente as instruções executadas por par e facilita o *pipelining*/vetorização ARM SVE no A64FX.
3. **Constante fora do ciclo e acumuladores locais:** remove ~25 milhões de multiplicações redundantes por `4 * epsilon` por passo e evita recargas de `r[i][0..2]` dentro do ciclo interno.

---

### 2.2. `computeAccelerations()` — Linhas 489–545 *(Segundo hotspot: ~33,1% na versão original)*

**O que mudou:**
- **Antes:** limpava `a[i][k]` num duplo ciclo `i, k`; percorria `i < N - 1`, `j = i + 1 .. N - 1`; calculava `rij[3]` e `rSqd` num ciclo `k = 0..2`; avaliava `f = 24 * (2 * pow(rSqd, -7) - pow(rSqd, -4))`; e num segundo ciclo `k = 0..2` atualizava `a[i][k] += rij[k] * f` e `a[j][k] -= rij[k] * f`.
- **Depois:**
  - Pré-calcula `const double c24 = 24.0;` fora dos ciclos e desenrola a inicialização `a[i][0..2] = 0.0`.
  - Mantém o ciclo triangular `i = 0 .. N - 2`, `j = i + 1 .. N - 1` (3.ª Lei de Newton, ótimo para execução sequencial).
  - Carrega `rxi, ryi, rzi` e os acumuladores `axi = a[i][0]`, `ayi = a[i][1]`, `azi = a[i][2]` em variáveis locais antes do ciclo em `j`, escrevendo de volta em `a[i][0..2]` apenas uma vez após terminar o ciclo em `j`.
  - Substitui `pow(rSqd, -4)` ($= r^{-8}$) e `pow(rSqd, -7)` ($= r^{-14}$) por multiplicações algébricas:
    ```cpp
    double inv_r2 = 1.0 / rSqd;
    double inv_r4 = inv_r2 * inv_r2;
    double inv_r8 = inv_r4 * inv_r4;
    double inv_r14 = inv_r8 * inv_r4 * inv_r2;
    double f = c24 * (2.0 * inv_r14 - inv_r8);
    ```
  - Calcula `fx = dx * f`, `fy = dy * f`, `fz = dz * f` uma única vez e aplica `axi += fx` e `a[j][0] -= fx` (e analogamente para `y`, `z`).

**Motivo físico e computacional:**
1. **Identidade algébrica exata:** como `rSqd` $= r^2$, tem-se $\text{rSqd}^{-4} = r^{-8} = \text{inv\_r8}$ e $\text{rSqd}^{-7} = r^{-14} = \text{inv\_r14}$. Eliminam-se $2 \times 12{,}5\text{ milhões} = 25\text{ milhões}$ de chamadas a `pow()` por invocação de `computeAccelerations()` (chamada 202 vezes numa corrida de Árgon).
2. **Redução de acessos à memória e multiplicações:** manter `rxi..rzi` e `axi..azi` em registos elimina $3$ leituras de `r[i]` e $3$ leituras + $3$ escritas de `a[i]` em cada iteração de `j`. Reutilizar `fx, fy, fz` evita calcular `rij[k] * f` duas vezes por componente (para `i` e para `j`).
3. **Nota de preparação para a Etapa 2 (OpenMP):** nesta 1.ª etapa sequencial, manter a 3.ª Lei de Newton ($j > i$) é a opção mais rápida porque executa metade das interações face à matriz completa $N \times N$. Na Etapa 2 (paralelização OpenMP), a atualização simultânea de `a[j][0..2]` introduz *data races*, pelo que será aí analisado o *trade-off* entre matriz completa $N \times N$, buffers privados por thread ou sincronização.

---

### 2.3. `VelocityVerlet()` — Linhas 548–593

**O que mudou:**
- **Antes:** tinha 3 ciclos separados sobre `i = 0 .. N - 1`, `j = 0 .. 2` (1.º meio-passo, 2.º meio-passo após `computeAccelerations()`, e verificação de paredes elásticas com dois `if` independentes) e recalculava `0.5 * a[i][j] * dt * dt`, `0.5 * a[i][j] * dt`, `2 * m * fabs(...) / dt` e `psum / (6 * L * L)`.
- **Depois:**
  - Pré-calcula `half_dt = 0.5 * dt`, `half_dt2 = 0.5 * dt * dt`, `wall_factor = 2.0 * m / dt` e `inv_area = 1.0 / (6.0 * L * L)` fora dos ciclos.
  - No 1.º ciclo, lê `aij = a[i][j]` uma só vez para atualizar `r[i][j]` e `v[i][j]`.
  - Após `computeAccelerations()`, **funde** a segunda atualização de velocidade com a verificação das paredes elásticas (`rij < 0.0 || rij >= L`) e a acumulação de `psum` num único ciclo, usando variáveis locais `vij` e `rij` e escrevendo `v[i][j]` apenas uma vez.
  - Devolve `psum * inv_area`.

**Motivo físico e computacional:**
1. **Fusão de ciclos (*loop fusion*):** reduz duas passagens completas sobre os arrays `v` e `r` a uma única passagem, melhorando a localidade temporal na cache L1 e reduzindo tráfego de memória.
2. **Condição mutuamente exclusiva:** como $L = V^{1/3} > 0$, uma coordenada não pode satisfazer simultaneamente $r_{ij} < 0$ e $r_{ij} \ge L$; combinar os dois testes em `rij < 0.0 || rij >= L` mantém o comportamento elástico idêntico evitando uma verificação redundante após colisão na parede esquerda.
3. **Constantes fora dos ciclos:** substitui divisões e multiplicações repetidas (`0.5 * dt`, `0.5 * dt * dt`, `2.0 * m / dt`, `1.0 / (6.0 * L * L)`) por constantes calculadas uma vez por passo.

---

### 2.4. `MeanSquaredVelocity()`, `Kinetic()`, `initializeVelocities()` e `main()`

- **`MeanSquaredVelocity()` (linhas 421–437):** pré-calcula `inv_N = 1.0 / N` e lê `v[i][0..2]` uma única vez em `vx, vy, vz`.
- **`Kinetic()` (linhas 440–452):** pré-calcula `half_m = 0.5 * m` fora do ciclo e acumula `v2_sum`, eliminando o ciclo interno `j = 0..2` e $N$ multiplicações por `m` e divisões por `2.0` em cada passo temporal.
- **`initializeVelocities()` (linhas 596–648):** pré-calcula `inv_Nm = 1.0 / (N * m)` e funde os ciclos de subtração da velocidade do centro de massa (`v[i][j] -= vCM[j]`) e cálculo da soma de quadrados (`vSqdSum`), poupando uma passagem completa sobre `v` (sem alterar a ordem das chamadas a `gaussdist()`/`rand()`, garantindo o mesmo estado inicial).
- **`main()` (linhas 282–348):** pré-calcula `VolSI`, `dt_SI`, `temp_scale`, `gc_scale` e `Z_scale` fora do ciclo temporal `i = 0 .. NumTime`, evitando recomputar expressões invariantes em cada iteração e no sumário final.

---

## 3. Checklist de Verificação e Testes no Deucalion (Antes da Etapa 2)

Antes de avançar para a paralelização OpenMP (Etapa 2), devem ser executadas as seguintes validações num nó ARM do Deucalion (`normal-arm`, `GCC/15.2.0`):

### 3.1. Compilação sequencial (sem OpenMP)
```bash
ml purge
ml GCC/15.2.0
gcc -O3 -ffast-math -march=native -g -fno-omit-frame-pointer src/MD1.cpp -lm -o md1
```
*(Opcionalmente, compilar também uma cópia do código original do professor para `md_orig` com exatamente o mesmo compilador/nó para comparação direta.)*

### 3.2. Paridade física e numérica dos outputs
Executar `./md1 < inputdata.txt` e validar:
1. **`cp_average.txt`:** confirmar que os valores médios termodinâmicos coincidem com a versão original do professor (tolerância esperada da ordem de $10^{-5}$ a $10^{-12}$ por reordenação de somas em vírgula flutuante, especialmente com `-ffast-math`):
   - Temperatura média: `~131.45506 K`
   - Pressão média: `~131001228.45076 Pa`
   - Compressibilidade: `Z ≈ 3.42449`
   - `N = 5000`, volume `V` e tempo total simulado idênticos.
2. **`cp_output.txt`:** confirmar que tem exatamente 202 linhas (cabeçalho + 201 passos `0..200`), que a evolução de `T(t)`, `P(t)`, `KE`, `PE` coincide com a baseline e que `Total En. (KE + PE)` se mantém fisicamente consistente/conservada (sem `NaN`/`Inf` nem divergência).
3. **`cp_traj.xyz`:** confirmar que mantém a mesma estrutura (`5000\n`, 6 bytes).

### 3.3. Medição de desempenho sequencial (1 thread) e novo profiling com `perf`
1. **Wall-clock time (1 core):**
   ```bash
   /usr/bin/time -p ./md1 < inputdata.txt
   ```
   - Comparar o tempo `real` obtido contra a baseline original (~`221.45 s`) e contra o valor de 1 thread do `MD.cpp` atual (`73.70 s` — nota: `MD1.cpp` deverá até ser **mais rápido** a 1 thread do que `73.70 s`, porque `computeAccelerations()` aqui usa o ciclo triangular $N(N-1)/2$ em vez da matriz completa $N \times N$!).
2. **Profiling pós-otimização sequencial:**
   ```bash
   perf record -g ./md1 < inputdata.txt
   perf report --stdio | head -n 40
   perf stat -e cycles,instructions,cache-references,cache-misses,branch-misses ./md1 < inputdata.txt
   ```
   - Verificar o novo peso percentual de `Potential()` vs `computeAccelerations()` (ambas agora com $N(N-1)/2$ pares e sem `pow`/`sqrt`), o ganho de IPC e a redução drástica de instruções/ciclos antes de desenhar a estratégia OpenMP da Etapa 2.
