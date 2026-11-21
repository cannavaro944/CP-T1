# Explicação Técnica — Etapa 2: Unificação das Otimizações Sequenciais (`MD1`) com OpenMP em `src/MD.cpp`

> **Âmbito desta etapa:** Unificação em `src/MD.cpp` das melhorias sequenciais desenvolvidas na Etapa 1 (`src/MD1.cpp`, documentadas em `EXPLICACAO_MD1_SEQUENCIAL.md`) com a paralelização OpenMP existente, garantindo **ausência de *data races***, preservação da física/saídas (`inputdata.txt`, `*_traj.xyz`, `*_output.txt`, `*_average.txt`) e correção de erros óbvios no `Makefile`.

---

## 1. Como as Melhorias de `MD1.cpp` e o OpenMP de `MD.cpp` Foram Unificados

| Função | Linhas em `src/MD.cpp` | Melhorias sequenciais herdadas de `MD1.cpp` | Paralelização OpenMP aplicada |
|---|---:|---|---|
| `main()` | 284–370 | Constantes invariantes pré-calculadas fora do ciclo temporal (`VolSI`, `dt_SI`, `temp_scale`, `gc_scale`, `Z_scale`). | **Sequencial** (dependência de estado entre passos temporais $t \to t+1$ e I/O ordenado). |
| `MeanSquaredVelocity()` | 443–457 | Constante `inv_N = 1.0 / N` fora do ciclo e carregamento único de `v[i][0..2]` em `vx, vy, vz`. | `#pragma omp parallel for reduction(+ : vx2, vy2, vz2) schedule(static)` |
| `Kinetic()` | 460–472 | Constante `half_m = 0.5 * m` fora do ciclo, acumulação direta em `v2_sum` (elimina ciclo `j = 0..2` e $N$ multiplicações/divisões por passo). | `#pragma omp parallel for reduction(+ : v2_sum) schedule(static)` |
| `Potential()` | 475–502 | Eliminação de `sqrt()`/`pow()` via `inv_r2` e `inv_r6`; ciclo triangular `i < N - 1`, `j = i + 1 .. N - 1`; constante `pot_factor = 8.0 * epsilon` (`2.0 * 4.0 * epsilon`); escalares `rxi..rzi` e `local_pot`. | `#pragma omp parallel for reduction(+ : Pot) schedule(dynamic, 16)` |
| `computeAccelerations()` | 507–547 | Substituição de `pow(rSqd, -7)` e `pow(rSqd, -4)` por `inv_r2`, `inv_r4`, `inv_r8`, `inv_r14`; constante `c24 = 24.0` fora do ciclo; `rxi..rzi` e acumuladores locais `ax, ay, az` em registos. | `#pragma omp parallel for schedule(static)` sobre a matriz completa $N \times N$ (sem escrita em `a[j]` e sem ciclo de limpeza prévio). |
| `VelocityVerlet()` | 550–585 | Constantes `half_dt`, `half_dt2`, `wall_factor = 2.0 * m / dt` e `inv_area = 1.0 / (6.0 * L * L)`; leitura única `aij`; fusão do 2.º meio-passo de velocidade com verificação de paredes (`rij < 0.0 \|\| rij >= L`) e `psum` usando escalares `vij, rij`. | 1.º ciclo: `#pragma omp parallel for schedule(static)`<br>2.º ciclo fundido: `#pragma omp parallel for reduction(+ : psum) schedule(static)` |
| `initialize()` / `initializeVelocities()` / `gaussdist()` | 399–440, 587–662 | Em `initializeVelocities()`: `inv_Nm = 1.0 / (N * m)` e fusão da subtração de `vCM` com o cálculo de `vSqdSum` numa só passagem sobre `v`. | **Sequencial** (executado uma única vez no arranque; `gaussdist()` usa `static bool available`, `static double gset` e `rand()`, que não são *thread-safe* e preservam o determinismo do estado inicial). |

---

## 2. Mapa das Diretivas OpenMP (`#pragma omp`, `reduction` e `schedule`)

Existem **6 regiões paralelas `#pragma omp parallel for`** em `src/MD.cpp`, todas dentro das funções invocadas a cada passo temporal:

1. **`MeanSquaredVelocity()` (linha 447):**
   ```cpp
   #pragma omp parallel for reduction(+ : vx2, vy2, vz2) schedule(static)
   ```
   - **Carga uniforme** ($O(N)$, 3 multiplicações e 3 somas por iteração `i`) $\Rightarrow$ `schedule(static)` minimiza o overhead de escalonamento.
   - `reduction(+ : vx2, vy2, vz2)` cria acumuladores privados por thread e soma-os de forma segura no fim da região.

2. **`Kinetic()` (linha 464):**
   ```cpp
   #pragma omp parallel for reduction(+ : v2_sum) schedule(static)
   ```
   - **Carga uniforme** ($O(N)$) $\Rightarrow$ `schedule(static)` com redução escalar única `v2_sum`, multiplicada uma vez por `half_m` após a barreira implícita.

3. **`Potential()` (linha 480):**
   ```cpp
   #pragma omp parallel for reduction(+ : Pot) schedule(dynamic, 16)
   ```
   - **Carga triangular desequilibrada** ($O(N^2 / 2)$): a iteração `i` executa `N - 1 - i` iterações internas em `j` (desde `4999` para `i = 0` até `1` para `i = 4998`).
   - Com `schedule(static)`, as primeiras threads receberiam muito mais trabalho do que as últimas. **`schedule(dynamic, 16)`** distribui blocos de 16 partículas dinamicamente, mantendo todas as threads ocupadas com baixo custo de sincronização.
   - Cada thread acumula primeiro em `local_pot` (registo) e adiciona à sua cópia privada da redução `Pot`, combinada no final da região.

4. **`computeAccelerations()` (linha 510):**
   ```cpp
   #pragma omp parallel for schedule(static)
   ```
   - **Carga perfeitamente uniforme** ($O(N^2)$): cada partícula `i` avalia `N - 1` interações (`j = 0 .. N - 1`, `j != i`).
   - **`schedule(static)`** atribui blocos contíguos e iguais de partículas `i` a cada thread, favorecendo localidade espacial na escrita final em `a[i][0..2]` e eliminando overhead de escalonamento dinâmico.

5. **`VelocityVerlet()` — 1.º ciclo (linha 558):**
   ```cpp
   #pragma omp parallel for schedule(static)
   ```
   - Atualiza `r[i][0..2]` e o 1.º meio-passo de `v[i][0..2]` com índices disjuntos por thread. A **barreira implícita** no final do ciclo garante que todas as posições `r` estão atualizadas antes de `computeAccelerations()` ser chamada.

6. **`VelocityVerlet()` — 2.º ciclo fundido (linha 570):**
   ```cpp
   #pragma omp parallel for reduction(+ : psum) schedule(static)
   ```
   - Atualiza o 2.º meio-passo de `v[i][0..2]`, testa colisão com as paredes elásticas (`rij < 0.0 || rij >= L`) e acumula a contribuição de momento em `psum` com `reduction(+ : psum)`.

---

## 3. Prevenção de *Data Races* e *Trade-off* em `computeAccelerations()`

### 3.1. O problema na versão sequencial (`MD1.cpp`)
Na Etapa 1 (`src/MD1.cpp`), `computeAccelerations()` tira partido da 3.ª Lei de Newton ($\mathbf{F}_{ji} = -\mathbf{F}_{ij}$) percorrendo apenas $i < j$:
```cpp
axi += fx;        // escreve no acumulador da partícula i
a[j][0] -= fx;    // ESCREVE na memória partilhada da partícula j!
```
Se este ciclo fosse paralelizado diretamente sobre `i`, múltiplas threads com diferentes `i` tentariam atualizar simultaneamente o mesmo `a[j][0..2]`, provocando ***data races*** graves (corrupção silenciosa das acelerações e divergência física) e forte *false sharing* nas linhas de cache de `a`.

### 3.2. Análise das alternativas em OpenMP e escolha aplicada
| Alternativa | Pares avaliados | Memória extra | Sincronização / Riscos | Decisão |
|---|---:|---|---|---|
| **1. Ciclo triangular + `#pragma omp atomic` em `a[j]`** | $N(N-1)/2$ | Nenhuma | Contenção atómica massiva ($3 \times 12{,}5\text{M}$ operações atómicas por passo) $\Rightarrow$ degrada severamente a escalabilidade em 48 cores. | Descartada |
| **2. Ciclo triangular + buffers privados `a_priv[tid][N][3]` + redução** | $N(N-1)/2$ | $p \times N \times 3 \times 8\text{ B}$ (~5,8 MB a 48T) | Sem data races no cálculo, mas exige limpeza e redução de 48 matrizes por passo + maior pressão de cache/memória HBM2 no A64FX e escalonamento dinâmico. | Alternativa viável, maior complexidade |
| **3. Matriz completa $N \times N$ por partícula `i` + acumuladores locais `ax, ay, az`** | $N(N-1)$ | Nenhuma | **Zero *data races***, **zero sincronização**, **zero ciclo de limpeza de `a`**, **zero *false sharing* no ciclo interno**, carga uniforme com `schedule(static)` e escrita única em `a[i][0..2]`. | **Adotada em `src/MD.cpp`** |

**Resumo do *trade-off*:** embora a abordagem por partícula `i` (matriz completa $N \times N$) dobre o número de pares avaliados em `computeAccelerations()` relativamente à 3.ª Lei de Newton de `MD1.cpp`, ela elimina completamente conflitos de escrita, dispensa o ciclo de inicialização de `a` a zero, mantém os acumuladores `ax, ay, az` em registos durante todo o ciclo interno de `j = 0 .. N - 1`, permite `schedule(static)` sem desbalanceamento e escala quase linearmente nos 48 cores ARM do Deucalion (no registo anterior `err_escala.e2206420`, atingiu **39,41×** de speedup e **1,87 s** a 48 threads).

### 3.3. Segurança em `psum`, `Pot`, `MeanSquaredVelocity` e `Kinetic`
- Em todas as reduções escalares (`psum`, `Pot`, `vx2/vy2/vz2`, `v2_sum`), o uso de `reduction(+ : ...)` garante cópias privadas inicializadas a `0.0` em cada thread e combinação determinística/segura na saída da região paralela.
- Todas as variáveis de ciclo (`int i`, `int j`) e auxiliares (`rxi..rzi`, `dx..dz`, `vij`, `rij`, `aij`, `local_pot`, `ax..az`) são declaradas dentro do corpo paralelo, sendo estritamente privadas a cada thread.
- Os arrays globais são lidos sem conflito (`r` é *read-only* em `Potential()` e `computeAccelerations()`) ou escritos em índices `i` estritamente disjuntos (`r[i]`, `v[i]`, `a[i]`).

---

## 4. Verificação e Correções no `Makefile` e `escalabilidade.sh`

- **`Makefile` (corrigido):**
  1. `clean:` alterado de `rm ./md` para `rm -f ./md`, eliminando o erro não-fatal (`rm: cannot remove './md': No such file or directory` / `Error 1`) que aparecia nas linhas 1–2 de `err_escala.e2206420`.
  2. `run:` alterado para `run: md`, garantindo que `make run` compila automaticamente `./md` antes de executar (requisito explícito do enunciado `WA1_26-27.pdf`).
  3. Dependência `$(SRC)/MD.cpp` normalizada para `$(SRC)MD.cpp` (evitando barra duplicada `src//MD.cpp`).
  4. Mantidas todas as flags (`-O3 -ffast-math -fopenmp -march=native -g -fno-omit-frame-pointer`) e `export OMP_NUM_THREADS ?= 48`.
- **`escalabilidade.sh` (verificado sem alterações):** já carrega `GCC/15.2.0`, corre `make clean && make` e testa `OMP_NUM_THREADS` em `{1, 2, 4, 8, 16, 32, 48}` com `/usr/bin/time -p ./md < inputdata.txt`. Com `rm -f ./md` no `Makefile`, passa a executar sem mensagens de erro iniciais em `err_escala.e%j`.

---

## 5. Pontos de Atenção para Testar no Deucalion

### 5.1. Compilação e teste rápido via `make run`
```bash
ml purge
ml GCC/15.2.0
make clean
make run
```
Confirmar que `make run` compila do zero sem warnings/erros e executa com 48 threads por omissão.

### 5.2. Paridade física e numérica (1 thread vs 48 threads vs `MD1.cpp`)
Executar com `OMP_NUM_THREADS=1` e `OMP_NUM_THREADS=48` e comparar `cp_average.txt` e `cp_output.txt`:
- Valores esperados de referência para Árgon (`N = 5000`, `100 K`, `35000 mol/m³`, 201 passos):
  - Temperatura média: `~131.45506 K`
  - Pressão média: `~131001228.45076 Pa`
  - Compressibilidade: `Z ≈ 3.42449`
- Verificar que a energia total (`KE + PE`) em `cp_output.txt` coincide com `MD1.cpp` (diferenças aceitáveis apenas ao nível de arredondamento de vírgula flutuante pela ordem de redução paralela e `-ffast-math`, sem desvios físicos nem `NaN`).

### 5.3. Escalabilidade forte (1 a 48 threads)
Submeter o script Slurm no nó ARM (`normal-arm`):
```bash
sbatch escalabilidade.sh
```
- Validar `resultados_tempos.txt` (`stdout`, com médias termodinâmicas por configuração) e `err_escala.e<jobid>` (`stderr`, com os 7 blocos `real/user/sys` de `/usr/bin/time -p`, agora sem erro do `make clean`).
- Recalcular `Speedup`, `Eficiência` e fração serial de **Karp-Flatt** ($e = \frac{1/S_p - 1/p}{1 - 1/p}$) e comparar com a corrida anterior (`73.70 s` em 1T $\to$ `1.87 s` em 48T).
- *(Opcional para o relatório)* Comparar também o tempo a **1 thread** de `MD1.cpp` (3.ª Lei de Newton, $N(N-1)/2$) com `MD.cpp` a 1 thread (matriz $N \times N$) para quantificar exatamente o custo inicial da duplicação de pares em `computeAccelerations()` face ao ganho obtido aos 48 cores.

### 5.4. Profiling com `perf` no Fujitsu A64FX
```bash
OMP_NUM_THREADS=48 perf record -g ./md < inputdata.txt
perf report --stdio | head -n 40
OMP_NUM_THREADS=48 perf stat -e cycles,instructions,cache-references,cache-misses,branch-misses ./md < inputdata.txt
```
Pontos a inspecionar no `perf`:
- Distribuição de tempo entre `computeAccelerations`, `Potential`, `VelocityVerlet` e runtime OpenMP (`libgomp` — tempo em barreiras/espera ativa).
- Possível impacto da arquitetura NUMA do A64FX (4 CMGs de 12 cores) ao passar de 32 para 48 threads; caso se queira testar afinidade de threads sem alterar código, experimentar `OMP_PROC_BIND=close OMP_PLACES=cores` vs configuração por omissão.
