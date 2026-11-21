# Contexto Técnico Completo — Computação Paralela (Fase 1: Dinâmica Molecular)

> **Objetivo deste ficheiro:** Servir de contexto técnico completo, autossuficiente e estruturado para revisão do relatório IEEE (2 páginas), preparação da defesa oral e explicação técnica do código de Dinâmica Molecular paralelizado com OpenMP.

---

=== VISAO_GERAL_E_REQUISITOS ===

## 1. Identificação e Ambiente de Execução
- **Unidade Curricular:** Computação Paralela (MEI — Universidade do Minho)
- **Autores / Grupo:** Pedro Lindo (`pg65436`) e Gonçalo Silva (`pg65412`)
- **Pasta de submissão esperada:** `pg65436_pg65412/` *(neste repositório Git local, os ficheiros estão diretamente na raiz `CP-T1/`)*
- **Problema:** Simulação de Dinâmica Molecular (*Molecular Dynamics — MD*) de um gás real de Árgon (`Ar`) com potencial de Lennard-Jones, integração numérica *Velocity Verlet* e condições de fronteira de paredes elásticas (baseado em *FoleyLab/MolecularDynamics*).
- **Hardware / Cluster:** Supercomputador **Deucalion**, partição `normal-arm` — nó ARM com processador **Fujitsu A64FX** (48 cores organizados em 4 CMGs de 12 cores, 32 GB HBM2, instruções vetoriais ARM SVE).
- **Compilador:** `GCC/15.2.0` (`ml purge && ml GCC/15.2.0`)
- **Flags atuais:** `-O3 -ffast-math -fopenmp -march=native -g -fno-omit-frame-pointer`
- **Profiling:** `perf` (`perf record` / `perf report` / `perf stat`)

## 2. Regras e Critérios de Avaliação do Enunciado (`extra/WA1_26-27.pdf`)
- **Prazo de entrega:** 25 de novembro, até às 23:59 (defesa oral na apresentação da Fase 2).
- **Formato da submissão:** ficheiro `.zip` que, ao descomprimir, gera a pasta `pg65436_pg65412/` contendo **apenas**:
  1. Relatório em **PDF** no template **IEEE**, com limite estrito de **2 páginas** (excluindo eventuais anexos).
  2. Subdiretoria `src/` com todo o código-fonte (sem executáveis nem ficheiros de saída/lixo).
  3. Um único `Makefile` na raiz que gera o executável `md` e disponibiliza o alvo `make run` para **compilar e executar** a versão paralela otimizada.
- **Critérios de avaliação:**
  - **30% — Tempo de execução** da implementação OpenMP nos nós ARM do Deucalion (via `make run`, podendo definir `OMP_NUM_THREADS ?= ...` no `Makefile`).
  - **30% — Abordagem de paralelização OpenMP e legibilidade do código** (*penalização forte se existirem data races*).
  - **40% — Qualidade do relatório** (análise de escalabilidade forte, profiling, modelos de desempenho e métricas que expliquem os resultados).

---

=== INVENTARIO_DA_PASTA ===

## 1. Árvore de Ficheiros do Repositório Atual (`CP-T1/`)

```text
CP-T1/
├── .DS_Store                       (6 148 B / 6.00 KiB)       [excluído: metadados macOS]
├── Makefile                        (279 B)                    [incluído na íntegra]
├── conference-template-a4.pdf      (226 428 B / 221.12 KiB)   [excluído: binário/PDF]
├── cp_average.txt                  (351 B)                    [excluído: dados brutos de saída]
├── cp_output.txt                   (25 052 B / 24.46 KiB)     [excluído: dados brutos de saída]
├── cp_traj.xyz                     (6 B)                      [excluído: dados brutos de saída]
├── err_escala.e2206420             (327 B)                    [incluído na íntegra]
├── escalabilidade.sh               (542 B)                    [incluído na íntegra]
├── inputdata.txt                   (22 B)                     [incluído na íntegra]
├── extra/
│   ├── Boas_Praticas.pdf           (12 792 967 B / 12.20 MiB) [excluído: binário/PDF]
│   ├── MD_simulations.pdf          (684 991 B / 668.94 KiB)   [excluído: binário/PDF]
│   ├── Profiling (1).pdf           (8 719 930 B / 8.32 MiB)   [excluído: binário/PDF]
│   ├── RELATORIO.docx              (24 355 B / 23.78 KiB)     [excluído: binário/DOCX — ver resumo]
│   └── WA1_26-27.pdf               (265 258 B / 259.04 KiB)   [excluído: binário/PDF — ver requisitos]
└── src/
    └── MD.cpp                      (21 111 B / 20.62 KiB)     [incluído na íntegra]
```

## 2. Para que serve cada ficheiro e pasta

### Código-fonte, Compilação e Execução (núcleo do projeto)
| Ficheiro / Pasta | Tamanho | Para que serve |
|---|---:|---|
| `src/MD.cpp` | 21 111 B | **Código-fonte principal** (653 linhas) do simulador de Dinâmica Molecular em C/C++ com paralelização OpenMP. Implementa a inicialização da rede cúbica e velocidades gaussianas, cálculo de forças de Lennard-Jones (`computeAccelerations`), energia potencial (`Potential`), integração temporal *Velocity Verlet* com paredes elásticas (`VelocityVerlet`) e médias termodinâmicas (`MeanSquaredVelocity`, `Kinetic`). |
| `Makefile` | 279 B | **Automatização de compilação e execução.** Compila `src/MD.cpp` com `gcc` e as flags de otimização/OpenMP/profiling para o binário `./md`, define `OMP_NUM_THREADS ?= 48`, disponibiliza `clean` e o alvo `run` (`./md < inputdata.txt`) usado na avaliação. |
| `inputdata.txt` | 22 B | **Ficheiro de entrada padrão (`stdin`)** do simulador. Passa os 4 parâmetros interativos lidos via `scanf`: prefixo dos ficheiros de saída (`cp`), gás nobre (`Ar`), temperatura inicial (`100` K) e densidade numérica (`35000` mol/m³). |
| `escalabilidade.sh` | 542 B | **Script Slurm para o Deucalion** (partição `normal-arm`, 1 nó, 48 CPUs, timeout 30 min). Carrega `GCC/15.2.0`, recompila o projeto e executa o teste de **escalabilidade forte** variando `OMP_NUM_THREADS` em `{1, 2, 4, 8, 16, 32, 48}` com `/usr/bin/time -p`. |

### Outputs e Logs gerados pelas execuções no Deucalion
| Ficheiro | Tamanho | Para que serve |
|---|---:|---|
| `err_escala.e2206420` | 327 B | **Ficheiro `stderr` do job Slurm `2206420`** (`#SBATCH --error=err_escala.e%j`). Contém um erro inicial não-fatal do `make clean` (porque `./md` ainda não existia) seguido dos **7 blocos de medição `real/user/sys`** produzidos por `/usr/bin/time -p` (que escreve para `stderr`). |
| `cp_average.txt` | 351 B | **Resumo termodinâmico final** gerado pelo simulador (`prefix + "_average.txt"`). Regista tempo total simulado, temperatura média `T`, pressão média `P`, `PV/nT`, fator de compressibilidade `Z`, volume `V` e `N`. Usado para **validação de correção física** (garantir que a versão paralela produz os mesmos resultados que a sequencial). |
| `cp_output.txt` | 25 052 B | **Série temporal instantânea** gerada pelo simulador (`prefix + "_output.txt"`). Contém cabeçalho + 201 linhas (passos `0..200`) com `time (s)`, `T(t)`, `P(t)`, energia cinética, energia potencial e energia total. Permite validar conservação de energia e estabilidade numérica passo a passo. |
| `cp_traj.xyz` | 6 B | **Ficheiro de trajetória para visualização (VMD)** (`prefix + "_traj.xyz"`). Tem apenas 6 bytes (`5000\n`) porque `main()` escreve o número de partículas no início, mas `VelocityVerlet()` recebe `FILE *fp` sem escrever coordenadas por passo (evitando I/O massivo). |

### Documentação, Template e Material de Apoio (não submeter no `.zip` final)
| Ficheiro / Pasta | Tamanho | Para que serve |
|---|---:|---|
| `conference-template-a4.pdf` | 226 428 B | **Template oficial IEEE (A4)** de referência visual/tipográfica para formatar o relatório final de 2 páginas. |
| `extra/WA1_26-27.pdf` | 265 258 B | **Enunciado oficial do trabalho (Work Assignment Phase 1)** com objetivos, metodologia recomendada, regras de submissão, anexos do `Makefile` e pesos da avaliação. |
| `extra/RELATORIO.docx` | 24 355 B | **Rascunho atual do relatório** em Word escrito pelo grupo (Pedro Lindo e Gonçalo Silva), com todas as secções (Abstract, Introdução, Perfilamento, Otimização/OpenMP, Escalabilidade Forte e Conclusão) pronto para revisão técnica e formatação final no template IEEE de 2 páginas. |
| `extra/MD_simulations.pdf` | 684 991 B | **Slides de apoio da UC** que explicam a física e o algoritmo base da simulação de Dinâmica Molecular (potencial de Lennard-Jones, Velocity Verlet e paredes elásticas). |
| `extra/Profiling (1).pdf` | 8 719 930 B | **Slides de apoio da UC** sobre metodologia de profiling (uso de `perf`, interpretação de hotspots, métricas de CPU/memória). |
| `extra/Boas_Praticas.pdf` | 12 792 967 B | **Material de apoio da UC** com boas práticas de programação paralela, análise de desempenho e escrita de relatórios técnicos. |
| `.DS_Store` | 6 148 B | Ficheiro oculto gerado automaticamente pelo macOS Finder (sem utilidade técnica; **deve ser removido** antes de criar o `.zip` de entrega). |

### Ficheiros ausentes na diretoria local
- **`resultados_tempos.txt`** — destino de `stdout` definido em `#SBATCH --output=resultados_tempos.txt` no `escalabilidade.sh` (conteria os separadores `=== Teste com $t threads ===` e o output de consola do `./md`, incluindo a confirmação das médias termodinâmicas em cada corrida). Não está presente localmente.
- **Outputs de `perf`** (`perf.data`, `perf report`, `perf stat`) — não estão guardados na pasta local (os números do profiling inicial foram transcritos para `extra/RELATORIO.docx`).
- **`projeto.zip` / pasta `pg65436_pg65412/`** — ainda não existem nesta diretoria; deverão ser criados no momento da submissão final.

---

=== ARQUITETURA_TECNICA_E_OTIMIZACOES ===

## 1. Fluxo de Execução e Complexidade de `src/MD.cpp`

Com `inputdata.txt` (`cp`, `Ar`, `100`, `35000`), a simulação corre com **`N = 5000` partículas** (fixo na linha 211) e **`NumTime = 200`** (total de **201 passos temporais**: `for (i = 0; i < NumTime + 1; i++)`).

```text
main()
 ├─ Lê stdin (prefix="cp", atype="Ar", Tinit=100 K, rho=35000 mol/m³)
 ├─ Calcula Vol, L = Vol^(1/3), dt = 0.5e-14 / timefac, abre cp_traj.xyz / cp_output.txt / cp_average.txt
 ├─ initialize()                 [Sequencial, O(N) — rede cúbica + initializeVelocities() com gaussdist()]
 ├─ computeAccelerations()       [Paralelo OpenMP, O(N²) — acelerações iniciais t=0]
 └─ Loop temporal i = 0..200     [Sequencial — dependência de estado entre passo t e t+1]
     ├─ VelocityVerlet(dt, ...)  [O(N) + chamada a computeAccelerations() O(N²) + O(N) c/ redução psum]
     ├─ MeanSquaredVelocity()    [Paralelo OpenMP, O(N) — redução vx2, vy2, vz2]
     ├─ Kinetic()                [Paralelo OpenMP, O(N) — redução kin]
     ├─ Potential()              [Paralelo OpenMP, O(N²/2) — ciclo triangular j > i, redução Pot]
     └─ fprintf(ofp, ...)        [Sequencial — I/O de 1 linha por passo + fflush(stdout)]
```

### Mapa das regiões OpenMP em `src/MD.cpp`
| Função | Linhas | Complexidade | Diretiva OpenMP | Estratégia / Justificação |
|---|---:|---|---|---|
| `Potential()` | 465–491 | $O(N^2 / 2)$ | `#pragma omp parallel for reduction(+ : Pot) schedule(dynamic, 16)` | Ciclo triangular (`j = i+1..N-1`) multiplicado por `2.0` no fim pela simetria dos pares; `schedule(dynamic, 16)` equilibra a carga decrescente por `i`; `reduction(+ : Pot)` evita data race. |
| `computeAccelerations()` | 496–536 | $O(N^2)$ | `#pragma omp parallel for schedule(static)` | Ciclo completo `i = 0..N-1`, `j = 0..N-1` (`if (i == j) continue;`). Abdica da simetria da 3.ª Lei de Newton para **eliminar data races** em `a[j][k]` sem atomics nem buffers privados. Carga uniforme permite `schedule(static)`. Acumuladores escalares `ax, ay, az` em registos eliminam o ciclo de limpeza de `a` e evitam *false sharing*. |
| `VelocityVerlet()` (passo 1) | 543–549 | $O(N)$ | `#pragma omp parallel for schedule(static)` | Atualiza posições `r[i][0..2]` e meio-passo de velocidades `v[i][0..2]` de forma independente por partícula `i`. |
| `VelocityVerlet()` (passo 2) | 554–568 | $O(N)$ | `#pragma omp parallel for reduction(+ : psum) schedule(static)` | Funde o segundo meio-passo de `v[i][j]` com a deteção/reflexão nas paredes elásticas (`r < 0` ou `r >= L`) e acumula momento transferido (`psum`) por redução. |
| `MeanSquaredVelocity()` | 436–446 | $O(N)$ | `#pragma omp parallel for reduction(+ : vx2, vy2, vz2) schedule(static)` | Soma dos quadrados das componentes de velocidade com redução paralela. |
| `Kinetic()` | 449–461 | $O(N)$ | `#pragma omp parallel for reduction(+ : kin) schedule(static)` | Soma da energia cinética de cada partícula com redução paralela. |

### Otimizações matemáticas e de compilação aplicadas
1. **Substituição de `pow()` e `sqrt()` por multiplicações diretas:**
   - Em `Potential()`: calcula `inv_r2 = 1.0 / rSqd` e `inv_r6 = inv_r2 * inv_r2 * inv_r2`, obtendo `inv_r6 * inv_r6 - inv_r6` ($r^{-12} - r^{-6}$).
   - Em `computeAccelerations()`: calcula `inv_r2`, `inv_r4`, `inv_r8`, `inv_r14 = inv_r8 * inv_r4 * inv_r2` e `f = 24.0 * (2.0 * inv_r14 - inv_r8)` (derivada do potencial de Lennard-Jones).
2. **Flags de compilação (`Makefile`):** `-O3 -ffast-math -fopenmp -march=native -g -fno-omit-frame-pointer` (ativam otimizações agressivas, reordenação de ponto flutuante, auto-vetorização para a microarquitetura ARM A64FX e símbolos/frame pointers para `perf`).
3. **Evolução registada no histórico Git:**
   - `8d56a89` (*Initial commit*): versão base do código com `Makefile` e `batchscript.sh`.
   - `9cd14d9` (*update md.cpp*): otimizações matemáticas + paralelização OpenMP.
   - `ea96d4f` (*excutavel*): substituição de `batchscript.sh` por `escalabilidade.sh`.
   - `36d1056` (*reverse*): revertida uma experiência anterior que usava `schedule(guided)` em `Potential()` (voltou a `schedule(dynamic, 16)`) e dividia o ciclo interno de `computeAccelerations()` em dois ciclos `0..i` e `i+1..N` (voltou a `0..N` com `if (i == j) continue;`).
   - `95fad0d` (*resultados*): gravação dos outputs `cp_*` e de `err_escala.e2206420`.

## 2. Dados Reportados no Rascunho do Relatório (`extra/RELATORIO.docx`)
*(Estes números constam do rascunho em Word e não aparecem nos ficheiros de log presentes na pasta — úteis para cruzamento na revisão):*
- **Profiling inicial (`perf record`, `N = 200`, nó A64FX):**
  - `Potential`: **65,4%** do tempo de CPU ($O(N^2)$)
  - `computeAccelerations`: **33,1%** do tempo de CPU ($O(N^2)$, ciclo triangular original)
  - Total dos dois hotspots: **> 98,5%**; rotinas $O(N)$ (`VelocityVerlet`, `MeanSquaredVelocity`, `Kinetic`): **< 1,5%**.
- **Evolução de desempenho (`N = 5000`, `201` passos):**
  - Baseline original (sequencial sem otimizações): **221,45 s**
  - Versão otimizada (1 thread, `-O3 -ffast-math -march=native` + simplificação algébrica): **73,70 s** (~3,00× mais rápida que a baseline)
  - Versão paralela preliminar referida no texto: **5,11 s**
  - Versão paralela final (48 threads): **1,87 s** (speedup paralelo **39,41×** vs 1T otimizado; redução total de **99,2%** vs baseline de 221,45 s).
- **Validação física/termodinâmica reportada (invariante entre 1 e 48 threads):**
  - Temperatura média: `131,45506 K`
  - Pressão média: `131001228,45076 Pa`
  - Fator de compressibilidade: `Z = 3,42449`

---

=== PONTOS_DE_ATENCAO_PARA_REVISAO ===

1. **Bug crítico no `Makefile` para a avaliação (`make run`):** o enunciado (`WA1_26-27.pdf`) diz explicitamente que `make run` deve *compilar e executar* o código. Atualmente o `Makefile` tem `run:` **sem dependência `md`** (deveria ser `run: md`). Se o avaliador executar `make run` numa pasta limpa, vai falhar porque `./md` ainda não foi compilado.
2. **Erro não-fatal no `Makefile` (`clean`):** a linha `rm ./md` falha com erro 1 se `./md` não existir (visível nas linhas 1–2 de `err_escala.e2206420`). Alterar para `rm -f ./md`.
3. **Pequena discrepância aritmética na Tabela I de `RELATORIO.docx`:** para `p = 4` threads, $T_1 / T_4 = 73.70 / 18.49 = 3.9859$, logo a eficiência é $3.9859 / 4 = \mathbf{99.65\%}$ (`99.6%` ou `99.7%`), e não `99.8%` (valor que resulta de dividir o speedup já arredondado `3.99 / 4`).
4. **Limpeza do `.zip` final de entrega:** garantir que o `.zip` gera a pasta `pg65436_pg65412/` contendo apenas o PDF IEEE (2 páginas), `Makefile` e `src/MD.cpp` — **excluir** `.DS_Store`, `extra/`, `cp_*`, `err_escala.*`, `escalabilidade.sh`, `contexto_projeto.md` e o binário `md`.
5. **Aspetos técnicos relevantes para reforçar no relatório (40% da nota):**
   - **Trade-off da matriz completa ($N \times N$) em `computeAccelerations()`:** calcula o dobro dos pares ($25\text{M}$ em vez de $12.5\text{M}$), mas ganha ausência total de sincronização, escrita única por partícula (`a[i][0..2]`), `schedule(static)` e potencial de auto-vetorização. Vale a pena justificar quantitativamente este trade-off (e mencionar se uma versão com acumulação privada ou divisão do loop `0..i` / `i+1..N` foi testada, dado o commit `36d1056`).
   - **Análise da queda de eficiência aos 48 cores (96,8% -> 82,1%):** além da Lei de Amdahl (fração sequencial + I/O + ~1000 barreiras implícitas OpenMP ao longo de 201 passos), notar a topologia do **Fujitsu A64FX** (4 CMGs NUMA de 12 cores cada), o aumento do tempo `user` de `73.68 s` (1T) para `81.70 s` (48T) indicando ~8 s de overhead paralelo/contenção, e a possibilidade de aplicar métricas formais (ex.: fração serial experimental de **Karp-Flatt**, $e = \frac{1/S_p - 1/p}{1 - 1/p}$).
   - **Nota sobre variáveis não utilizadas / vestígios em `MD.cpp`:** `double F[MAXPART][3]` (linha 60) e `char trash[10000]` / `FILE *infp` (linhas 93–94) são vestígios do código original; `MeanSquaredVelocity()` e `Kinetic()` percorrem ambas `v[i][0..2]` separadamente (poderiam até partilhar o cálculo, embora sejam $< 1.5\%$ do tempo).

---

=== Makefile ===

```makefile
CC = gcc
SRC = src/
CFLAGS = -O3 -ffast-math -fopenmp -march=native -g -fno-omit-frame-pointer

.DEFAULT_GOAL = md

md:     $(SRC)/MD.cpp
	$(CC) $(CFLAGS) $(SRC)MD.cpp -lm -o md    

clean:
	rm ./md    

export OMP_NUM_THREADS ?= 48
run:
	./md     < inputdata.txt
```

---

=== escalabilidade.sh ===

```bash
#!/bin/bash
#SBATCH -A f202500010hpcvlabuminhoa
#SBATCH -p normal-arm
#SBATCH -t 30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --output=resultados_tempos.txt
#SBATCH --error=err_escala.e%j

# Preparar ambiente
ml purge
ml GCC/15.2.0
make clean
make

# Testar com diferentes números de threads e medir o tempo (wall-clock time)
for t in 1 2 4 8 16 32 48; do
    export OMP_NUM_THREADS=$t
    echo "=== Teste com $t threads ==="
    /usr/bin/time -p ./md < inputdata.txt
    echo ""
done
```

---

=== inputdata.txt ===

*(4 linhas úteis / 5 linhas totais incluindo newline final; incluído na íntegra)*

```text
cp
Ar
100
35000
```

- Linha 1 (`cp`): prefixo dos ficheiros de saída (`cp_traj.xyz`, `cp_output.txt`, `cp_average.txt`).
- Linha 2 (`Ar`): gás nobre simulado (Árgon; seleciona fatores de conversão SI e `NumTime = 200`, `dt = 0.5e-14 / timefac`).
- Linha 3 (`100`): temperatura inicial $T_{\text{init}} = 100\text{ K}$.
- Linha 4 (`35000`): densidade numérica $\rho = 35000\text{ mol/m}^3$.

---

=== src/MD.cpp ===

*(653 linhas; incluído na íntegra)*

```cpp
/*
 MD.c - a simple molecular dynamics program for simulating real gas properties
 of Lennard-Jones particles.

 Copyright (C) 2016  Jonathan J. Foley IV, Chelsea Sweet, Oyewumi Akinfenwa

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 Electronic Contact:  foleyj10@wpunj.edu
 Mail Contact:   Prof. Jonathan Foley
 Department of Chemistry, William Paterson University
 300 Pompton Road
 Wayne NJ 07470

 */
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Number of particles
int N;

//  Lennard-Jones parameters in natural units!
double sigma = 1.;
double epsilon = 1.;
double m = 1.;
double kB = 1.;

double NA = 6.022140857e23;
double kBSI = 1.38064852e-23; // m^2*kg/(s^2*K)

//  Size of box, which will be specified in natural units
double L;

//  Initial Temperature in Natural Units
double Tinit; // 2;
//  Vectors!
//
const int MAXPART = 5001;
//  Position
double r[MAXPART][3];
//  Velocity
double v[MAXPART][3];
//  Acceleration
double a[MAXPART][3];
//  Force
double F[MAXPART][3];

// atom type
char atype[10];
//  Function prototypes
//  initialize positions on simple cubic lattice, also calls function to
//  initialize velocities
void initialize();
//  update positions and velocities using Velocity Verlet algorithm
//  print particle coordinates to file for rendering via VMD or other animation
//  software return 'instantaneous pressure'
double VelocityVerlet(double dt, int iter, FILE *fp);
//  Compute Force using F = -dV/dr
//  solve F = ma for use in Velocity Verlet
void computeAccelerations();
//  Numerical Recipes function for generation gaussian distribution
double gaussdist();
//  Initialize velocities according to user-supplied initial Temperature (Tinit)
void initializeVelocities();
//  Compute total potential energy from particle coordinates
double Potential();
//  Compute mean squared velocity from particle velocities
double MeanSquaredVelocity();
//  Compute total kinetic energy from particle mass and velocities
double Kinetic();

int main() {

  //  variable delcarations
  int i;
  double dt, Vol, Temp, Press, Pavg, Tavg, rho;
  double VolFac, TempFac, PressFac, timefac;
  double KE, PE, mvs, gc, Z;
  char trash[10000], prefix[1000], tfn[1000], ofn[1000], afn[1000];
  FILE *infp, *tfp, *ofp, *afp;

  printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  printf("                  WELCOME TO WILLY P CHEM MD!\n");
  printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  printf("\n  ENTER A TITLE FOR YOUR CALCULATION!\n");
  scanf("%s", prefix);
  strcpy(tfn, prefix);
  strcat(tfn, "_traj.xyz");
  strcpy(ofn, prefix);
  strcat(ofn, "_output.txt");
  strcpy(afn, prefix);
  strcat(afn, "_average.txt");

  printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  printf("                  TITLE ENTERED AS '%s'\n", prefix);
  printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");

  /*     Table of values for Argon relating natural units to SI units:
   *     These are derived from Lennard-Jones parameters from the article
   *     "Liquid argon: Monte carlo and molecular dynamics calculations"
   *     J.A. Barker , R.A. Fisher & R.O. Watts
   *     Mol. Phys., Vol. 21, 657-673 (1971)
   *
   *     mass:     6.633e-26 kg          = one natural unit of mass for argon,
   *by definition energy:   1.96183e-21 J      = one natural unit of energy for
   *argon, directly from L-J parameters length:   3.3605e-10  m         = one
   *natural unit of length for argon, directly from L-J parameters
   *     volume:   3.79499-29 m^3        = one natural unit of volume for argon,
   *by length^3 time:     1.951e-12 s           = one natural unit of time for
   *argon, by length*sqrt(mass/energy)
   ***************************************************************************************/

  //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //  Edit these factors to be computed in terms of basic properties in natural
  //  units of the gas being simulated

  printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  printf("  WHICH NOBLE GAS WOULD YOU LIKE TO SIMULATE? (DEFAULT IS ARGON)\n");
  printf("\n  FOR HELIUM,  TYPE 'He' THEN PRESS 'return' TO CONTINUE\n");
  printf("  FOR NEON,    TYPE 'Ne' THEN PRESS 'return' TO CONTINUE\n");
  printf("  FOR ARGON,   TYPE 'Ar' THEN PRESS 'return' TO CONTINUE\n");
  printf("  FOR KRYPTON, TYPE 'Kr' THEN PRESS 'return' TO CONTINUE\n");
  printf("  FOR XENON,   TYPE 'Xe' THEN PRESS 'return' TO CONTINUE\n");
  printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  scanf("%s", atype);

  if (strcmp(atype, "He") == 0) {

    VolFac = 1.8399744000000005e-29;
    PressFac = 8152287.336171632;
    TempFac = 10.864459551225972;
    timefac = 1.7572698825166272e-12;

  } else if (strcmp(atype, "Ne") == 0) {

    VolFac = 2.0570823999999997e-29;
    PressFac = 27223022.27659913;
    TempFac = 40.560648991243625;
    timefac = 2.1192341945685407e-12;

  } else if (strcmp(atype, "Ar") == 0) {

    VolFac = 3.7949992920124995e-29;
    PressFac = 51695201.06691862;
    TempFac = 142.0950000000000;
    timefac = 2.09618e-12;
    // strcpy(atype,"Ar");

  } else if (strcmp(atype, "Kr") == 0) {

    VolFac = 4.5882712000000004e-29;
    PressFac = 59935428.40275003;
    TempFac = 199.1817584391428;
    timefac = 8.051563913585078e-13;

  } else if (strcmp(atype, "Xe") == 0) {

    VolFac = 5.4872e-29;
    PressFac = 70527773.72794868;
    TempFac = 280.30305642163006;
    timefac = 9.018957925790732e-13;

  } else {

    VolFac = 3.7949992920124995e-29;
    PressFac = 51695201.06691862;
    TempFac = 142.0950000000000;
    timefac = 2.09618e-12;
    strcpy(atype, "Ar");
  }
  printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  printf("\n                     YOU ARE SIMULATING %s GAS! \n", atype);
  printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");

  printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  printf("\n  YOU WILL NOW ENTER A FEW SIMULATION PARAMETERS\n");
  printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  printf("\n\n  ENTER THE INTIAL TEMPERATURE OF YOUR GAS IN KELVIN\n");
  scanf("%lf", &Tinit);
  // Make sure temperature is a positive number!
  if (Tinit < 0.) {
    printf("\n  !!!!! ABSOLUTE TEMPERATURE MUST BE A POSITIVE NUMBER!  PLEASE "
           "TRY AGAIN WITH A POSITIVE TEMPERATURE!!!\n");
    exit(0);
  }
  // Convert initial temperature from kelvin to natural units
  Tinit /= TempFac;

  printf("\n\n  ENTER THE NUMBER DENSITY IN moles/m^3\n");
  printf("  FOR REFERENCE, NUMBER DENSITY OF AN IDEAL GAS AT STP IS ABOUT 40 "
         "moles/m^3\n");
  printf("  NUMBER DENSITY OF LIQUID ARGON AT 1 ATM AND 87 K IS ABOUT 35000 "
         "moles/m^3\n");

  scanf("%lf", &rho);

  N = 5000;
  Vol = N / (rho * NA);

  Vol /= VolFac;

  //  Limiting N to MAXPART for practical reasons
  if (N >= MAXPART) {

    printf("\n\n\n  MAXIMUM NUMBER OF PARTICLES IS %i\n\n  PLEASE ADJUST YOUR "
           "INPUT FILE ACCORDINGLY \n\n",
           MAXPART);
    exit(0);
  }
  //  Check to see if the volume makes sense - is it too small?
  //  Remember VDW radius of the particles is 1 natural unit of length
  //  and volume = L*L*L, so if V = N*L*L*L = N, then all the particles
  //  will be initialized with an interparticle separation equal to 2xVDW radius
  if (Vol < N) {

    printf("\n\n\n  YOUR DENSITY IS VERY HIGH!\n\n");
    printf("  THE NUMBER OF PARTICLES IS %i AND THE AVAILABLE VOLUME IS %f "
           "NATURAL UNITS\n",
           N, Vol);
    printf("  SIMULATIONS WITH DENSITY GREATER THAN 1 PARTCICLE/(1 Natural "
           "Unit of Volume) MAY DIVERGE\n");
    printf("  PLEASE ADJUST YOUR INPUT FILE ACCORDINGLY AND RETRY\n\n");
    exit(0);
  }
  // Vol = L*L*L;
  // Length of the box in natural units:
  L = pow(Vol, (1. / 3));

  //  Files that we can write different quantities to
  tfp = fopen(tfn, "w"); //  The MD trajectory, coordinates of every particle at
                         //  each timestep
  ofp = fopen(
      ofn,
      "w"); //  Output of other quantities (T, P, gc, etc) at every timestep
  afp = fopen(afn, "w"); //  Average T, P, gc, etc from the simulation

  int NumTime;
  if (strcmp(atype, "He") == 0) {

    // dt in natural units of time s.t. in SI it is 5 f.s. for all other gasses
    dt = 0.2e-14 / timefac;
    //  We will run the simulation for NumTime timesteps.
    //  The total time will be NumTime*dt in natural units
    //  And NumTime*dt multiplied by the appropriate conversion factor for time
    //  in seconds
    NumTime = 50000;
  } else {
    dt = 0.5e-14 / timefac;
    NumTime = 200;
  }

  //  Put all the atoms in simple crystal lattice and give them random
  //  velocities that corresponds to the initial temperature we have specified
  initialize();

  //  Based on their positions, calculate the ininial intermolecular forces
  //  The accellerations of each particle will be defined from the forces and
  //  their mass, and this will allow us to update their positions via Newton's
  //  law
  computeAccelerations();

  // Print number of particles to the trajectory file
  fprintf(tfp, "%i\n", N);

  //  We want to calculate the average Temperature and Pressure for the
  //  simulation The variables need to be set to zero initially
  Pavg = 0;
  Tavg = 0;

  int tenp = floor(NumTime / 10);
  fprintf(ofp,
          "  time (s)              T(t) (K)              P(t) (Pa)           "
          "Kinetic En. (n.u.)     Potential En. (n.u.) Total En. (n.u.)\n");
  printf("  PERCENTAGE OF CALCULATION COMPLETE:\n  [");
  for (i = 0; i < NumTime + 1; i++) {

    //  This just prints updates on progress of the calculation for the users
    //  convenience
    if (i == tenp)
      printf(" 10 |");
    else if (i == 2 * tenp)
      printf(" 20 |");
    else if (i == 3 * tenp)
      printf(" 30 |");
    else if (i == 4 * tenp)
      printf(" 40 |");
    else if (i == 5 * tenp)
      printf(" 50 |");
    else if (i == 6 * tenp)
      printf(" 60 |");
    else if (i == 7 * tenp)
      printf(" 70 |");
    else if (i == 8 * tenp)
      printf(" 80 |");
    else if (i == 9 * tenp)
      printf(" 90 |");
    else if (i == 10 * tenp)
      printf(" 100 ]\n");
    fflush(stdout);

    // This updates the positions and velocities using Newton's Laws
    // Also computes the Pressure as the sum of momentum changes from wall
    // collisions / timestep which is a Kinetic Theory of gasses concept of
    // Pressure
    Press = VelocityVerlet(dt, i + 1, tfp);
    Press *= PressFac;

    //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //  Now we would like to calculate somethings about the system:
    //  Instantaneous mean velocity squared, Temperature, Pressure
    //  Potential, and Kinetic Energy
    //  We would also like to use the IGL to try to see if we can extract the
    //  gas constant
    mvs = MeanSquaredVelocity();
    KE = Kinetic();
    PE = Potential();

    // Temperature from Kinetic Theory
    Temp = m * mvs / (3 * kB) * TempFac;

    // Instantaneous gas constant and compressibility - not well defined because
    // pressure may be zero in some instances because there will be zero wall
    // collisions, pressure may be very high in some instances because there
    // will be a number of collisions
    gc = NA * Press * (Vol * VolFac) / (N * Temp);
    Z = Press * (Vol * VolFac) / (N * kBSI * Temp);

    Tavg += Temp;
    Pavg += Press;

    fprintf(ofp, "  %8.4e  %20.8f  %20.8f %20.8f  %20.8f  %20.8f \n",
            i * dt * timefac, Temp, Press, KE, PE, KE + PE);
  }

  // Because we have calculated the instantaneous temperature and pressure,
  // we can take the average over the whole simulation here
  Pavg /= NumTime;
  Tavg /= NumTime;
  Z = Pavg * (Vol * VolFac) / (N * kBSI * Tavg);
  gc = NA * Pavg * (Vol * VolFac) / (N * Tavg);
  fprintf(afp, "  Total Time (s)      T (K)               P (Pa)      PV/nT "
               "(J/(mol K))         Z           V (m^3)              N\n");
  fprintf(afp,
          " --------------   -----------        ---------------   "
          "--------------   ---------------   ------------   -----------\n");
  fprintf(afp,
          "  %8.4e  %15.5f       %15.5f     %10.5f       %10.5f        %10.5e  "
          "       %i\n",
          i * dt * timefac, Tavg, Pavg, gc, Z, Vol * VolFac, N);

  printf("\n  TO ANIMATE YOUR SIMULATION, OPEN THE FILE \n  '%s' WITH VMD "
         "AFTER THE SIMULATION COMPLETES\n",
         tfn);
  printf(
      "\n  TO ANALYZE INSTANTANEOUS DATA ABOUT YOUR MOLECULE, OPEN THE FILE \n "
      " '%s' WITH YOUR FAVORITE TEXT EDITOR OR IMPORT THE DATA INTO EXCEL\n",
      ofn);
  printf("\n  THE FOLLOWING THERMODYNAMIC AVERAGES WILL BE COMPUTED AND "
         "WRITTEN TO THE FILE  \n  '%s':\n",
         afn);
  printf("\n  AVERAGE TEMPERATURE (K):                 %15.5f\n", Tavg);
  printf("\n  AVERAGE PRESSURE  (Pa):                  %15.5f\n", Pavg);
  printf("\n  PV/nT (J * mol^-1 K^-1):                 %15.5f\n", gc);
  printf("\n  PERCENT ERROR of pV/nT AND GAS CONSTANT: %15.5f\n",
         100 * fabs(gc - 8.3144598) / 8.3144598);
  printf("\n  THE COMPRESSIBILITY (unitless):          %15.5f \n", Z);
  printf("\n  TOTAL VOLUME (m^3):                      %10.5e \n",
         Vol * VolFac);
  printf("\n  NUMBER OF PARTICLES (unitless):          %i \n", N);

  fclose(tfp);
  fclose(ofp);
  fclose(afp);

  return 0;
}

void initialize() {
  int n, p, i, j, k;
  double pos;

  // Number of atoms in each direction
  n = int(ceil(pow(N, 1.0 / 3)));

  //  spacing between atoms along a given direction
  pos = L / n;

  //  index for number of particles assigned positions
  p = 0;
  //  initialize positions
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      for (k = 0; k < n; k++) {
        if (p < N) {

          r[p][0] = (i + 0.5) * pos;
          r[p][1] = (j + 0.5) * pos;
          r[p][2] = (k + 0.5) * pos;
        }
        p++;
      }
    }
  }

  // Call function to initialize velocities
  initializeVelocities();

  /***********************************************
   *   Uncomment if you want to see what the initial positions and velocities
   are printf("  Printing initial positions!\n"); for (i=0; i<N; i++) { printf("
   %6.3e  %6.3e  %6.3e\n",r[i][0],r[i][1],r[i][2]);
   }

   printf("  Printing initial velocities!\n");
   for (i=0; i<N; i++) {
   printf("  %6.3e  %6.3e  %6.3e\n",v[i][0],v[i][1],v[i][2]);
   }
   */
}

//  Function to calculate the averaged velocity squared
double MeanSquaredVelocity() {
  double vx2 = 0, vy2 = 0, vz2 = 0;

#pragma omp parallel for reduction(+ : vx2, vy2, vz2) schedule(static)
  for (int i = 0; i < N; i++) {
    vx2 += v[i][0] * v[i][0];
    vy2 += v[i][1] * v[i][1];
    vz2 += v[i][2] * v[i][2];
  }
  return (vx2 + vy2 + vz2) / N;
}

//  Function to calculate the kinetic energy of the system
double Kinetic() {
  double kin = 0.;

#pragma omp parallel for reduction(+ : kin) schedule(static)
  for (int i = 0; i < N; i++) {
    double v2 = 0.;
    for (int j = 0; j < 3; j++) {
      v2 += v[i][j] * v[i][j];
    }
    kin += m * v2 / 2.;
  }
  return kin;
}

// Function to calculate the potential energy of the system

double Potential() {
  double Pot = 0.0;
  int i;

#pragma omp parallel for reduction(+ : Pot) schedule(dynamic, 16)
  for (i = 0; i < N - 1; i++) {
    double rxi = r[i][0];
    double ryi = r[i][1];
    double rzi = r[i][2];
    double local_pot = 0.0;

    for (int j = i + 1; j < N; j++) {
      double dx = rxi - r[j][0];
      double dy = ryi - r[j][1];
      double dz = rzi - r[j][2];

      double rSqd = dx * dx + dy * dy + dz * dz;
      double inv_r2 = 1.0 / rSqd;
      double inv_r6 = inv_r2 * inv_r2 * inv_r2;

      local_pot += (inv_r6 * inv_r6 - inv_r6);
    }
    Pot += local_pot;
  }
  // Multiplica por 2.0 devido à simetria dos pares (i,j) e (j,i)
  return 2.0 * 4.0 * epsilon * Pot;
}

//   Uses the derivative of the Lennard-Jones potential to calculate
//   the forces on each atom.  Then uses a = F/m to calculate the
//   accelleration of each atom.
void computeAccelerations() {
  int i;

#pragma omp parallel for schedule(static)
  for (i = 0; i < N; i++) {
    double rxi = r[i][0];
    double ryi = r[i][1];
    double rzi = r[i][2];

    // Acumuladores locais em registo (evita False Sharing na matriz 'a' e
    // dispensa o ciclo de limpeza)
    double ax = 0.0, ay = 0.0, az = 0.0;

    for (int j = 0; j < N; j++) {
      if (i == j)
        continue;

      double dx = rxi - r[j][0];
      double dy = ryi - r[j][1];
      double dz = rzi - r[j][2];

      double rSqd = dx * dx + dy * dy + dz * dz;

      double inv_r2 = 1.0 / rSqd;
      double inv_r4 = inv_r2 * inv_r2;
      double inv_r8 = inv_r4 * inv_r4;
      double inv_r14 = inv_r8 * inv_r4 * inv_r2;

      double f = 24.0 * (2.0 * inv_r14 - inv_r8);

      ax += dx * f;
      ay += dy * f;
      az += dz * f;
    }

    // Escreve na memória partilhada apenas uma vez no final
    a[i][0] = ax;
    a[i][1] = ay;
    a[i][2] = az;
  }
}

// returns sum of dv/dt*m/A (aka Pressure) from elastic collisions with walls
double VelocityVerlet(double dt, int iter, FILE *fp) {
  int i, j;
  double psum = 0.;

#pragma omp parallel for schedule(static)
  for (i = 0; i < N; i++) {
    for (j = 0; j < 3; j++) {
      r[i][j] += v[i][j] * dt + 0.5 * a[i][j] * dt * dt;
      v[i][j] += 0.5 * a[i][j] * dt;
    }
  }

  computeAccelerations();

// Funde a atualização de velocidade com a verificação das paredes elásticas
#pragma omp parallel for reduction(+ : psum) schedule(static)
  for (i = 0; i < N; i++) {
    for (j = 0; j < 3; j++) {
      v[i][j] += 0.5 * a[i][j] * dt;

      if (r[i][j] < 0.) {
        v[i][j] *= -1.;
        psum += 2 * m * fabs(v[i][j]) / dt;
      }
      if (r[i][j] >= L) {
        v[i][j] *= -1.;
        psum += 2 * m * fabs(v[i][j]) / dt;
      }
    }
  }

  return psum / (6 * L * L);
}

void initializeVelocities() {

  int i, j;

  for (i = 0; i < N; i++) {

    for (j = 0; j < 3; j++) {
      //  Pull a number from a Gaussian Distribution
      v[i][j] = gaussdist();
    }
  }

  // Vcm = sum_i^N  m*v_i/  sum_i^N  M
  // Compute center-of-mas velocity according to the formula above
  double vCM[3] = {0, 0, 0};

  for (i = 0; i < N; i++) {
    for (j = 0; j < 3; j++) {

      vCM[j] += m * v[i][j];
    }
  }

  for (i = 0; i < 3; i++)
    vCM[i] /= N * m;

  //  Subtract out the center-of-mass velocity from the
  //  velocity of each particle... effectively set the
  //  center of mass velocity to zero so that the system does
  //  not drift in space!
  for (i = 0; i < N; i++) {
    for (j = 0; j < 3; j++) {

      v[i][j] -= vCM[j];
    }
  }

  //  Now we want to scale the average velocity of the system
  //  by a factor which is consistent with our initial temperature, Tinit
  double vSqdSum, lambda;
  vSqdSum = 0.;
  for (i = 0; i < N; i++) {
    for (j = 0; j < 3; j++) {

      vSqdSum += v[i][j] * v[i][j];
    }
  }

  lambda = sqrt(3 * (N - 1) * Tinit / vSqdSum);

  for (i = 0; i < N; i++) {
    for (j = 0; j < 3; j++) {

      v[i][j] *= lambda;
    }
  }
}

//  Numerical recipes Gaussian distribution number generator
double gaussdist() {
  static bool available = false;
  static double gset;
  double fac, rsq, v1, v2;
  if (!available) {
    do {
      v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
      v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
      rsq = v1 * v1 + v2 * v2;
    } while (rsq >= 1.0 || rsq == 0.0);

    fac = sqrt(-2.0 * log(rsq) / rsq);
    gset = v1 * fac;
    available = true;

    return v2 * fac;
  } else {

    available = false;
    return gset;
  }
}
```

---

=== resultados_tempos.txt ===

```text
[FICHEIRO NÃO ENCONTRADO NA DIRETORIA ATUAL]
```

> **Nota:** `escalabilidade.sh` direciona `stdout` para `resultados_tempos.txt` (`#SBATCH --output=resultados_tempos.txt`), onde ficariam os cabeçalhos `=== Teste com $t threads ===` e o texto impresso pelo programa `./md`. Como `/usr/bin/time -p` escreve as medições para `stderr`, todos os tempos (`real`, `user`, `sys`) foram capturados em `err_escala.e2206420` (secção seguinte).

---

=== err_escala.e2206420 ===

```text
rm: cannot remove './md': No such file or directory
make: *** [Makefile:11: clean] Error 1
real 73.70
user 73.68
sys 0.00
real 37.02
user 73.85
sys 0.10
real 18.49
user 73.88
sys 0.01
real 9.26
user 73.93
sys 0.03
real 4.65
user 74.17
sys 0.00
real 2.38
user 75.29
sys 0.01
real 1.87
user 81.70
sys 0.12
```

### Tabela consolidada de escalabilidade forte (`err_escala.e2206420`, $N = 5000$, 201 passos)

| Threads ($p$) | `real` (s) | `user` (s) | `sys` (s) | Speedup ($S_p = T_1 / T_p$) | Eficiência ($E_p = S_p / p$) | Fração serial Karp-Flatt ($e$) |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 73.70 | 73.68 | 0.00 | 1.00× | 100.00% | — |
| 2 | 37.02 | 73.85 | 0.10 | 1.99× | 99.54% | 0.0046 |
| 4 | 18.49 | 73.88 | 0.01 | 3.99× | 99.65% | 0.0012 |
| 8 | 9.26 | 73.93 | 0.03 | 7.96× | 99.49% | 0.0007 |
| 16 | 4.65 | 74.17 | 0.00 | 15.85× | 99.06% | 0.0006 |
| 32 | 2.38 | 75.29 | 0.01 | 30.97× | 96.77% | 0.0011 |
| 48 | 1.87 | 81.70 | 0.12 | 39.41× | 82.11% | 0.0046 |

*(Fórmula de Karp-Flatt: $e = \frac{1/S_p - 1/p}{1 - 1/p}$; o salto de $e$ de `0.0011` em 32T para `0.0046` em 48T — acompanhado pelo salto de tempo `user` de `75.29 s` para `81.70 s` — evidencia aumento de overhead paralelo/contenção de memória/sincronização ao ocupar os 48 cores do A64FX, além da fração puramente sequencial prevista pela Lei de Amdahl.)*
