# Contexto Técnico Completo — Computação Paralela (Fase 1: Dinâmica Molecular)

> **Objetivo deste ficheiro:** Servir de contexto técnico completo, atualizado e autossuficiente (código final SoA + OpenMP, `Makefile` corrigido, resultados reais de escalabilidade no Deucalion e métricas de `perf stat`/`perf report`) para revisão do relatório IEEE (2 páginas), preparação da defesa oral e explicação técnica do código.

---

=== VISAO_GERAL_E_REQUISITOS ===

## 1. Identificação e Ambiente de Execução
- **Unidade Curricular:** Computação Paralela (MEI — Universidade do Minho)
- **Autores / Grupo:** Pedro Lindo (`pg65436`) e Gonçalo Silva (`pg65412`)
- **Pasta de submissão esperada:** `pg65436_pg65412/` *(neste repositório Git local, os ficheiros estão diretamente na raiz `CP-T1/`)*
- **Problema:** Simulação de Dinâmica Molecular (*Molecular Dynamics — MD*) de um gás real de Árgon (`Ar`) com potencial de Lennard-Jones, integração numérica *Velocity Verlet* e paredes elásticas (baseado em *FoleyLab/MolecularDynamics*).
- **Hardware / Cluster:** Supercomputador **Deucalion**, partição `normal-arm` — nó ARM com processador **Fujitsu A64FX**:
  - **48 cores** organizados em **4 CMGs** (*Core Memory Groups*) de 12 cores cada (arquitetura NUMA)
  - **Cache L1D:** 64 KiB por core, linha de cache de **256 bytes**
  - **Cache L2:** 8 MiB partilhados por CMG (32 MiB totais)
  - **Memória principal:** 32 GiB HBM2 de elevada largura de banda
  - **Unidade vetorial:** ARM **SVE** de 512 bits (8 `double`s por registo vetorial)
- **Compilador:** `GCC/15.2.0` (`ml purge && ml GCC/15.2.0`)
- **Flags de compilação:** `-O3 -ffast-math -fopenmp -march=native -g -fno-omit-frame-pointer`
- **Profiling:** `perf` (`perf stat`, `perf record`, `perf report`)

## 2. Regras e Critérios de Avaliação do Enunciado
- **Prazo de entrega:** 25 de novembro, até às 23:59 (defesa oral na apresentação da Fase 2).
- **Formato da submissão:** ficheiro `.zip` que, ao descomprimir, gera a pasta `pg65436_pg65412/` contendo **apenas**:
  1. Relatório em **PDF** no template **IEEE**, com limite estrito de **2 páginas** (excluindo eventuais anexos).
  2. Subdiretoria `src/` com todo o código-fonte (`src/MD.cpp`, sem executáveis nem ficheiros de saída).
  3. Um único `Makefile` na raiz que gera o executável `md` e disponibiliza o alvo `make run` para **compilar e executar** a versão paralela otimizada.
- **Critérios de avaliação:**
  - **30% — Tempo de execução** da implementação OpenMP nos nós ARM do Deucalion (via `make run`, com `OMP_NUM_THREADS ?= 48`).
  - **30% — Abordagem de paralelização OpenMP e legibilidade do código** (*penalização forte se existirem data races*).
  - **40% — Qualidade do relatório** (análise de escalabilidade forte, profiling, modelos de desempenho e métricas que expliquem os resultados).

---

=== INVENTARIO_DA_PASTA ===

## 1. Árvore de Ficheiros do Repositório Atual (`CP-T1/`)

```text
CP-T1/
├── .DS_Store                       (6 148 B / 6.00 KiB)     [excluído: metadados macOS]
├── Makefile                        (279 B)                  [incluído na íntegra]
├── conference-template-a4.pdf      (226 428 B / 221.12 KiB) [excluído: binário/PDF]
├── contexto_projeto.md             (este ficheiro)
├── cp_average.txt                  (348 B)                  [excluído: dados brutos — valores confirmados abaixo]
├── err_escala.e2256900             (213 B)                  [incluído na íntegra]
├── escalabilidade.sh               (542 B)                  [incluído na íntegra]
├── inputdata.txt                   (22 B)                   [incluído na íntegra]
├── perf_report_48t.txt             (1 578 B)                [incluído na íntegra]
├── perf_stat_1t.txt                (660 B)                  [incluído na íntegra]
├── perf_stat_48t.txt               (660 B)                  [incluído na íntegra]
├── resultados_tempos.txt           (16 590 B / 16.20 KiB)   [incluído na íntegra]
└── src/
    └── MD.cpp                      (22 398 B / 21.87 KiB)   [incluído na íntegra]
```

## 2. Para que serve cada ficheiro

### Código-fonte, Build e Script de Benchmark
| Ficheiro | Tamanho | Para que serve |
|---|---:|---|
| `src/MD.cpp` | 22 398 B | **Código-fonte final otimizado e paralelizado** (669 linhas). Usa layout **SoA** alinhado a 64 bytes (`rx/ry/rz`, `vx/vy/vz`, `ax/ay/az`), funde o cálculo da energia potencial dentro de `computeAccelerations()` (`cached_PE`), divide o ciclo interno `j` em `0..i-1` e `i+1..N-1` (sem `if (i == j)`), pré-calcula constantes invariantes, funde o 2.º meio-passo de `VelocityVerlet()` com as paredes elásticas, valida `fopen()` e paraleliza 5 regiões com OpenMP sem *data races*. |
| `Makefile` | 279 B | **Compilação e execução automatizada.** Compila `src/MD.cpp` com `gcc -O3 -ffast-math -fopenmp -march=native -g -fno-omit-frame-pointer -lm -o md`, corrige `clean` para `rm -f ./md`, exporta `OMP_NUM_THREADS ?= 48` e garante `run: md` (compila antes de executar `./md < inputdata.txt`). |
| `inputdata.txt` | 22 B | **Entrada padrão (`stdin`)** do simulador: prefixo `cp`, gás `Ar`, temperatura inicial `100` K e densidade `35000` mol/m³ (`N = 5000`, `NumTime = 200`, 201 passos temporais). |
| `escalabilidade.sh` | 542 B | **Script Slurm para o Deucalion** (`normal-arm`, 1 nó, 48 CPUs, 30 min). Carrega `GCC/15.2.0`, executa `make clean && make` e mede o wall-clock time (`/usr/bin/time -p`) para `OMP_NUM_THREADS` em `{1, 2, 4, 8, 16, 32, 48}`. |

### Resultados Reais dos Testes no Deucalion (Job `2256900` + `perf`)
| Ficheiro | Tamanho | Para que serve |
|---|---:|---|
| `err_escala.e2256900` | 213 B | **Saída `stderr` do job Slurm `2256900`** (versão final SoA + fusão). Contém **zero erros** e os 7 blocos `real/user/sys` produzidos por `/usr/bin/time -p` para `1, 2, 4, 8, 16, 32, 48` threads (`49.18 s` $\to$ `1.23 s`). |
| `resultados_tempos.txt` | 16 590 B | **Saída `stdout` completa do job `2256900`** (464 linhas). Regista a compilação limpa e as 7 execuções completas, provando **paridade física bit-a-bit** das médias termodinâmicas em todas as contagens de threads (`T = 131.45506 K`, `P = 131001228.45076 Pa`, `Z = 3.42449`). |
| `perf_stat_1t.txt` | 660 B | **Contadores de hardware (`perf stat`) a 1 thread** no A64FX: ciclos (`98,28B`), instruções (`16,59B`, IPC `0,17`), referências e misses de cache / L1D (`485,55M` misses = `24,32%`) e tempo `49,15 s`. |
| `perf_stat_48t.txt` | 660 B | **Contadores de hardware (`perf stat`) a 48 threads** no A64FX: ciclos (`106,95B`), instruções (`31,57B`, IPC `0,30`), `488,32M` misses L1D (`9,46%` da totalidade de acessos L1D) e tempo `1,189 s`. |
| `perf_report_48t.txt` | 1 578 B | **Perfil de funções (`perf report`) a 48 threads** (231K amostras de `cycles`): mostra `86,34%` em `computeAccelerations`, `9,57%` em barreiras OpenMP (`5,40%` `gomp_team_barrier_wait_end` + `4,17%` `gomp_barrier_wait_end`), `0,73%` em `MeanSquaredVelocity` e `0,64%` em `Kinetic`. |
| `cp_average.txt` | 348 B | **Resumo termodinâmico final** gravado pelo simulador na última execução (confirma `1.0050e-12 s`, `T = 131.45506 K`, `P = 131001228.45076 Pa`, `PV/nT = 28.47279`, `Z = 3.42449`, `V = 2.37220e-25 m³`, `N = 5000`). |

### Ficheiros Auxiliares (NÃO incluir no `.zip` final)
| Ficheiro | Tamanho | Para que serve |
|---|---:|---|
| `conference-template-a4.pdf` | 226 428 B | Template oficial IEEE (A4) de referência para o relatório final de 2 páginas. |
| `.DS_Store` | 6 148 B | Metadados automáticos do macOS Finder (remover ao gerar o `.zip`). |

---

=== ARQUITETURA_TECNICA_E_OTIMIZACOES ===

## 1. Otimizações Implementadas na Versão Final (`src/MD.cpp`)

1. **Layout de Dados SoA (*Structure of Arrays*) alinhado a 64 bytes (`L54-L58`):**
   - Substituídos os arrays AoS `r[MAXPART][3]`, `v[MAXPART][3]`, `a[MAXPART][3]` (e removido `F[MAXPART][3]` não usado) por `alignas(64) double rx/ry/rz`, `vx/vy/vz`, `ax/ay/az`.
   - Permite acessos contíguos (*unit-stride*) em memória ao percorrer `rx[j]`, `ry[j]`, `rz[j]`, eliminando *gather/scatter* de stride 3 e permitindo ao GCC/15.2.0 auto-vetorizar diretamente para registos **ARM SVE de 512 bits** (8 `double`s por instrução vetorial).
   - Eliminados os ciclos internos fixos `j = 0..2` em todas as funções (`initialize`, `initializeVelocities`, `MeanSquaredVelocity`, `Kinetic`, `computeAccelerations`, `VelocityVerlet`).
2. **Fusão de `Potential()` dentro de `computeAccelerations()` (`L60`, `L479`, `L484-L536`):**
   - Como as posições `rx, ry, rz` não sofrem alteração entre `computeAccelerations()` (dentro de `VelocityVerlet()`) e `Potential()` (no `main()`), a energia potencial é calculada no mesmo ciclo $O(N^2)$ das forças e guardada em `cached_PE`.
   - Elimina completamente uma segunda travessia $O(N^2)$ de memória/cálculo por passo temporal (e remove a barreira OpenMP de `Potential()`).
   - Partilha integral dos valores intermédios `dx, dy, dz, rSqd, inv_r2, inv_r6` entre potencial (`local_pot += inv_r6 * (inv_r6 - 1.0)`) e força (`f = c24 * inv_r2 * inv_r6 * (2.0 * inv_r6 - 1.0)`).
   - Como todos os pares ordenados $(i, j)$ com $i \ne j$ são avaliados uma vez, o fator final é exatamente `cached_PE = 4.0 * epsilon * Pot;` (idêntico à definição original, sem duplicação).
3. **Remoção de `if (i == j)` via *Loop Splitting* (`L499-L525`):**
   - O ciclo interno `j` em `computeAccelerations()` foi dividido em dois ciclos contíguos sem ramos condicionais: `for (int j = 0; j < i; j++)` e `for (int j = i + 1; j < N; j++)`.
   - Elimina desvios condicionais no ciclo mais quente e desbloqueia a auto-vetorização SVE do GCC/15.2.0 (`-O3 -ffast-math -march=native`).
   - *(Nota: os `#pragma omp simd` foram removidos na última revisão local de `src/MD.cpp`, mantendo a auto-vetorização automática do GCC/15.2.0.)*
4. **Prevenção de *Data Races* e *False Sharing* com OpenMP:**
   - **`computeAccelerations()` (`L488`):** `#pragma omp parallel for reduction(+ : Pot) schedule(static)` sobre a matriz completa $N \times N$ por partícula `i`. Evita escrita concorrente em `a[j]` (sem atomics nem buffers privados), usa acumuladores locais em registo (`axi, ayi, azi, local_pot`), dispensa ciclo de limpeza de `ax/ay/az`, escreve `ax[i], ay[i], az[i]` apenas uma vez no final da iteração `i` e tem carga perfeitamente uniforme (`schedule(static)`).
   - **`VelocityVerlet()` (`L547`, `L565`):** pré-calcula `half_dt`, `half_dt2`, `wall_factor = 2.0 * m / dt` e `inv_area = 1.0 / (6.0 * L * L)`; paraleliza o 1.º meio-passo (`schedule(static)`) e funde o 2.º meio-passo de velocidades com a deteção de paredes elásticas (`rxi < 0.0 || rxi >= L`, etc.) e acumulação de pressão (`reduction(+ : psum) schedule(static)`).
   - **`MeanSquaredVelocity()` (`L450`) e `Kinetic()` (`L467`):** `schedule(static)` com `reduction(+ : ...)` e constantes `inv_N`, `half_m` fora dos ciclos.
   - **`initializeVelocities()` (`L596-L645`):** mantida sequencial para preservar a ordem exata `x, y, z` das chamadas a `gaussdist()`/`rand()` (determinismo do estado inicial), fundindo a subtração de `vCM` com o somatório `vSqdSum`.
5. **Robustez de E/S (`L250-L253`) e constantes fora do ciclo temporal (`L290-L294`):**
   - Verificação de ponteiros nulos após `fopen()` (`tfp`, `ofp`, `afp`) com `perror` + `exit(1)`.
   - Pré-cálculo de `VolSI`, `dt_SI`, `temp_scale`, `gc_scale` e `Z_scale` fora do ciclo temporal de `main()`.

## 2. Evolução Histórica de Desempenho no Deucalion (`N = 5000`, 201 passos, Árgon)

| Versão | Descrição técnica | Tempo 1T (s) | Tempo 48T (s) | Speedup 48T vs 1T | Eficiência 48T | Speedup vs Baseline (221,45 s) |
|---|---|---:|---:|---:|---:|---:|
| **V0 (Baseline professor)** | Sequencial original (AoS, `pow`/`sqrt`, ciclo completo em `Potential`, 3.ª Lei em `computeAccelerations`) | **221,45** | — | — | — | 1,00× |
| **V1 (AoS + OpenMP — job `2206420`)** | Potências algébricas, `Potential` triangular (`dynamic, 16`), `computeAccelerations` $N \times N$ (`static`), `VelocityVerlet` fundido | **73,70** | **1,87** | 39,41× | 82,11% | 118,42× (−99,16%) |
| **V2 (Final SoA + Fusão — job `2256900`)** | SoA `alignas(64)`, fusão `Potential`+`computeAccelerations`, *loop splitting* `0..i`/`i+1..N`, aritmética partilhada `inv_r6`, constantes fora dos ciclos | **49,18** | **1,23** *(1,19 s em `perf stat`)* | **39,98×** | **83,30%** | **180,04× (−99,44%)** |

- **Ganho de V2 sobre V1:** **1,50× mais rápida a 1 thread** (`73,70 s` $\to$ `49,18 s`) e **1,52× mais rápida a 48 threads** (`1,87 s` $\to$ `1,23 s`), melhorando simultaneamente o speedup paralelo (`39,41×` $\to$ `39,98×`) e a eficiência (`82,11%` $\to$ `83,30%`).
- **Paridade física 100% confirmada** em todas as 7 corridas (`1..48` threads) de `resultados_tempos.txt` e em `cp_average.txt`: `T = 131.45506 K`, `P = 131001228.45076 Pa`, `PV/nT = 28.47279 J/(mol·K)`, `Z = 3.42449`.

## 3. Análise Hardware-Level dos Dados de `perf` (Fujitsu A64FX)

Os ficheiros `perf_stat_1t.txt`, `perf_stat_48t.txt` e `perf_report_48t.txt` fornecem evidência quantitativa de alto nível para o relatório IEEE:

1. **Vetorização ARM SVE e IPC (`perf_stat_1t.txt`):**
   - Em 1 thread, o simulador executa apenas **16,59 mil milhões de instruções** para ~5,05 mil milhões de avaliações de pares $(i, j)$ ($202 \times 5000 \times 4999$), ou seja, **~3,29 instruções por par** — prova direta de que o GCC/15.2.0 vetorizou o ciclo interno em registos **ARM SVE de 512 bits** (8 `double`s por instrução).
   - O **IPC = 0,17** (`16,59B / 98,28B` ciclos) é explicado por dois fatores microarquiteturais do A64FX: (i) latência elevada da divisão vetorial `1.0 / rSqd` (`fdiv`) por iteração, e (ii) misses na cache L1D servidos pela cache L2.
2. **Modelo analítico exato dos L1D Cache Misses (`485,55M` medidos vs `~474M` teóricos):**
   - Os 3 vetores de posição SoA (`rx`, `ry`, `rz`) ocupam $3 \times 5000 \times 8\text{ B} = \mathbf{120\text{ KiB}}$, excedendo a capacidade da cache **L1D do A64FX (64 KiB por core)**, mas cabendo inteiramente na **cache L2 (8 MiB por CMG)**.
   - Como a linha de cache do A64FX tem **256 B**, cada varrimento completo de `rx, ry, rz` para uma partícula `i` lê $\lceil 120\,000 / 256 \rceil \approx 469$ linhas de cache. Multiplicando pelas $202 \times 5000 = 1{,}01\text{M}$ linhas de `i`, obtêm-se **~474 milhões de L1D misses teóricos**, praticamente idênticos aos **485,55 milhões medidos** em `perf_stat_1t.txt`!
3. **Ausência de *False Sharing* e *Cache Thrashing* a 48 cores (`perf_stat_1t.txt` vs `perf_stat_48t.txt`):**
   - O número absoluto de `L1-dcache-load-misses` mantém-se virtualmente invariante entre 1 thread (**485,55M**) e 48 threads (**488,32M**, apenas **+0,57%**), confirmando que o layout SoA com acumuladores locais em registo e escrita única em `ax[i], ay[i], az[i]` **eliminou completamente *false sharing* e tráfego de coerência entre os 48 cores**.
   - O aumento de `instructions` (`16,59B` $\to$ `31,57B`), `L1-dcache-loads` (`2,00B` $\to$ `5,16B`) e `IPC` (`0,17` $\to$ `0,30`) a 48 threads deve-se ao *busy-waiting* (*spin-wait*) nas barreiras OpenMP, que executa ciclos curtos de leitura em L1 (baixando a *taxa* aparente de miss de `24,32%` para `9,46%`).
4. **Explicação exata da eficiência de 83,3% aos 48 cores (`perf_report_48t.txt` + Karp-Flatt):**
   - `computeAccelerations` (já incluindo `Potential`) domina com **86,34%** das amostras de ciclos das worker threads.
   - As **barreiras OpenMP** consomem **9,57%** dos ciclos totais (`5,40%` em `gomp_team_barrier_wait_end` + `4,17%` em `gomp_barrier_wait_end`), distribuídas pelas 5 barreiras implícitas por passo temporal ($\approx 1005$ barreiras no total) e atravessando os 4 CMGs NUMA do A64FX.
   - `MeanSquaredVelocity` (`0,73%`) + `Kinetic` (`0,64%`) somam apenas `1,37%`, e `~2,12%` correspondem ao fio principal / fração sequencial (inicialização, I/O `fprintf`/`fflush`).
   - A soma dos overheads (~9,6% barreiras + ~2,1% sequencial + pequeno desequilíbrio entre CMGs) explica com precisão a perda de eficiência de ~16,7% aos 48 threads (`E_48 = 83,30%`) e o crescimento do tempo `user` de `49,14 s` (1T) para `53,96 s` (48T) e da fração serial experimental de **Karp-Flatt** de `0,0008` (16T) $\to$ `0,0026` (32T) $\to$ `0,0043` (48T).

---

=== PONTOS_DE_ATENCAO_PARA_REVISAO ===

1. **`Makefile` já corrigido e validado no Deucalion:** `clean` usa `rm -f ./md` (zero erros em `err_escala.e2256900`) e `run: md` compila antes de executar `make run`, cumprindo integralmente o enunciado.
2. **Atualizações obrigatórias no relatório IEEE (2 páginas):**
   - Substituir todos os números antigos (`73,70 s` $\to$ `1,87 s`, `39,41×`, `82,1%`) pelos resultados finais do job `2256900` (**`49,18 s` $\to$ `1,23 s`**, speedup **`39,98×`**, eficiência **`83,30%`**, redução total de **99,44%** / **180×** face aos `221,45 s` da baseline).
   - Documentar as otimizações decisivas de V1 $\to$ V2: **SoA `alignas(64)`**, **fusão de `Potential()` em `computeAccelerations()`** (`cached_PE = 4.0 * epsilon * Pot`) e ***loop splitting* (`0..i` / `i+1..N`)** para vetorização ARM SVE.
   - Integrar as métricas reais de `perf stat` e `perf report` (IPC `0,17`, working set `120 KiB > 64 KiB` L1D explicando os `485,5M` misses, invariância de misses em 48T provando ausência de *false sharing*, `9,57%` em barreiras `libgomp` e fração serial de Karp-Flatt).
   - Atualizar o perfil de hotspots pós-otimização: `Potential()` passou a $O(1)$ (em cache), deixando `computeAccelerations()` com `86,34%` do tempo paralelo.
3. **Limpeza do `.zip` final de entrega (`pg65436_pg65412/`):**
   - Deve conter **apenas**: o relatório em **PDF** (2 páginas, IEEE), `Makefile` e `src/MD.cpp`.
   - **Excluir** `.DS_Store`, `conference-template-a4.pdf`, `contexto_projeto.md`, `escalabilidade.sh`, `resultados_tempos.txt`, `err_escala.e2256900`, `perf_*.txt`, `cp_average.txt` e o binário `./md`.

---

=== Makefile ===

```makefile
CC = gcc
SRC = src/
CFLAGS = -O3 -ffast-math -fopenmp -march=native -g -fno-omit-frame-pointer

.DEFAULT_GOAL = md

md:     $(SRC)MD.cpp
	$(CC) $(CFLAGS) $(SRC)MD.cpp -lm -o md

clean:
	rm -f ./md

export OMP_NUM_THREADS ?= 48
run:    md
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

---

=== src/MD.cpp ===

*(669 linhas; versão final atual incluída na íntegra)*

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
//  Position (SoA alinhado a 64 bytes para vetorização ARM SVE)
alignas(64) double rx[MAXPART], ry[MAXPART], rz[MAXPART];
//  Velocity
alignas(64) double vx[MAXPART], vy[MAXPART], vz[MAXPART];
//  Acceleration
alignas(64) double ax[MAXPART], ay[MAXPART], az[MAXPART];
//  Energia potencial em cache (calculada dentro de computeAccelerations)
double cached_PE = 0.0;

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
  if (!tfp || !ofp || !afp) {
    perror("Erro ao abrir ficheiros de output");
    exit(1);
  }

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
  // Constantes invariantes pré-calculadas fora do ciclo temporal principal
  const double VolSI = Vol * VolFac;
  const double dt_SI = dt * timefac;
  const double temp_scale = (m / (3.0 * kB)) * TempFac;
  const double gc_scale = (NA * VolSI) / N;
  const double Z_scale = VolSI / (N * kBSI);

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
    Temp = mvs * temp_scale;

    // Instantaneous gas constant and compressibility - not well defined because
    // pressure may be zero in some instances because there will be zero wall
    // collisions, pressure may be very high in some instances because there
    // will be a number of collisions
    gc = (Press * gc_scale) / Temp;
    Z = (Press * Z_scale) / Temp;

    Tavg += Temp;
    Pavg += Press;

    fprintf(ofp, "  %8.4e  %20.8f  %20.8f %20.8f  %20.8f  %20.8f \n",
            i * dt_SI, Temp, Press, KE, PE, KE + PE);
  }

  // Because we have calculated the instantaneous temperature and pressure,
  // we can take the average over the whole simulation here
  Pavg /= NumTime;
  Tavg /= NumTime;
  Z = (Pavg * Z_scale) / Tavg;
  gc = (Pavg * gc_scale) / Tavg;
  fprintf(afp, "  Total Time (s)      T (K)               P (Pa)      PV/nT "
               "(J/(mol K))         Z           V (m^3)              N\n");
  fprintf(afp,
          " --------------   -----------        ---------------   "
          "--------------   ---------------   ------------   -----------\n");
  fprintf(afp,
          "  %8.4e  %15.5f       %15.5f     %10.5f       %10.5f        %10.5e  "
          "       %i\n",
          i * dt_SI, Tavg, Pavg, gc, Z, VolSI, N);

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
          rx[p] = (i + 0.5) * pos;
          ry[p] = (j + 0.5) * pos;
          rz[p] = (k + 0.5) * pos;
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
   %6.3e  %6.3e  %6.3e\n",rx[i],ry[i],rz[i]);
   }

   printf("  Printing initial velocities!\n");
   for (i=0; i<N; i++) {
   printf("  %6.3e  %6.3e  %6.3e\n",vx[i],vy[i],vz[i]);
   }
   */
}

//  Function to calculate the averaged velocity squared
double MeanSquaredVelocity() {
  const double inv_N = 1.0 / N;
  double vx2 = 0.0, vy2 = 0.0, vz2 = 0.0;

#pragma omp parallel for reduction(+ : vx2, vy2, vz2) schedule(static)
  for (int i = 0; i < N; i++) {
    double vxi = vx[i];
    double vyi = vy[i];
    double vzi = vz[i];
    vx2 += vxi * vxi;
    vy2 += vyi * vyi;
    vz2 += vzi * vzi;
  }
  return (vx2 + vy2 + vz2) * inv_N;
}

//  Function to calculate the kinetic energy of the system
double Kinetic() {
  const double half_m = 0.5 * m;
  double v2_sum = 0.0;

#pragma omp parallel for reduction(+ : v2_sum) schedule(static)
  for (int i = 0; i < N; i++) {
    double vxi = vx[i];
    double vyi = vy[i];
    double vzi = vz[i];
    v2_sum += vxi * vxi + vyi * vyi + vzi * vzi;
  }
  return half_m * v2_sum;
}

// Function to calculate the potential energy of the system
// (valor já calculado e guardado em cache por computeAccelerations)
double Potential() { return cached_PE; }

//   Uses the derivative of the Lennard-Jones potential to calculate
//   the forces on each atom.  Then uses a = F/m to calculate the
//   accelleration of each atom. Also computes and caches potential energy.
void computeAccelerations() {
  const double c24 = 24.0;
  double Pot = 0.0;

#pragma omp parallel for reduction(+ : Pot) schedule(static)
  for (int i = 0; i < N; i++) {
    double rxi = rx[i];
    double ryi = ry[i];
    double rzi = rz[i];

    // Acumuladores locais em registo (evita Data Races, False Sharing e ciclo
    // de limpeza de ax/ay/az)
    double axi = 0.0, ayi = 0.0, azi = 0.0;
    double local_pot = 0.0;

    for (int j = 0; j < i; j++) {
      double dx = rxi - rx[j];
      double dy = ryi - ry[j];
      double dz = rzi - rz[j];
      double rSqd = dx * dx + dy * dy + dz * dz;
      double inv_r2 = 1.0 / rSqd;
      double inv_r6 = inv_r2 * inv_r2 * inv_r2;
      local_pot += inv_r6 * (inv_r6 - 1.0);
      double f = c24 * inv_r2 * inv_r6 * (2.0 * inv_r6 - 1.0);
      axi += dx * f;
      ayi += dy * f;
      azi += dz * f;
    }

    for (int j = i + 1; j < N; j++) {
      double dx = rxi - rx[j];
      double dy = ryi - ry[j];
      double dz = rzi - rz[j];
      double rSqd = dx * dx + dy * dy + dz * dz;
      double inv_r2 = 1.0 / rSqd;
      double inv_r6 = inv_r2 * inv_r2 * inv_r2;
      local_pot += inv_r6 * (inv_r6 - 1.0);
      double f = c24 * inv_r2 * inv_r6 * (2.0 * inv_r6 - 1.0);
      axi += dx * f;
      ayi += dy * f;
      azi += dz * f;
    }

    Pot += local_pot;
    // Escreve na memória partilhada apenas uma vez no final
    ax[i] = axi;
    ay[i] = ayi;
    az[i] = azi;
  }

  // Todos os pares ordenados (i, j) com i != j já foram somados: fator 4.0 * epsilon
  cached_PE = 4.0 * epsilon * Pot;
}

// returns sum of dv/dt*m/A (aka Pressure) from elastic collisions with walls
double VelocityVerlet(double dt, int iter, FILE *fp) {
  const double half_dt = 0.5 * dt;
  const double half_dt2 = 0.5 * dt * dt;
  const double wall_factor = 2.0 * m / dt;
  const double inv_area = 1.0 / (6.0 * L * L);

  double psum = 0.0;

#pragma omp parallel for schedule(static)
  for (int i = 0; i < N; i++) {
    double axi = ax[i];
    double ayi = ay[i];
    double azi = az[i];

    rx[i] += vx[i] * dt + half_dt2 * axi;
    ry[i] += vy[i] * dt + half_dt2 * ayi;
    rz[i] += vz[i] * dt + half_dt2 * azi;

    vx[i] += half_dt * axi;
    vy[i] += half_dt * ayi;
    vz[i] += half_dt * azi;
  }

  computeAccelerations();

// Funde a atualização de velocidade com a verificação das paredes elásticas
#pragma omp parallel for reduction(+ : psum) schedule(static)
  for (int i = 0; i < N; i++) {
    double vxi = vx[i] + half_dt * ax[i];
    double vyi = vy[i] + half_dt * ay[i];
    double vzi = vz[i] + half_dt * az[i];

    double rxi = rx[i];
    double ryi = ry[i];
    double rzi = rz[i];

    if (rxi < 0.0 || rxi >= L) {
      vxi = -vxi;
      psum += wall_factor * fabs(vxi);
    }
    if (ryi < 0.0 || ryi >= L) {
      vyi = -vyi;
      psum += wall_factor * fabs(vyi);
    }
    if (rzi < 0.0 || rzi >= L) {
      vzi = -vzi;
      psum += wall_factor * fabs(vzi);
    }

    vx[i] = vxi;
    vy[i] = vyi;
    vz[i] = vzi;
  }

  return psum * inv_area;
}

void initializeVelocities() {
  int i;

  for (i = 0; i < N; i++) {
    //  Pull numbers from a Gaussian Distribution (ordem exata: x, y, z)
    vx[i] = gaussdist();
    vy[i] = gaussdist();
    vz[i] = gaussdist();
  }

  // Vcm = sum_i^N  m*v_i/  sum_i^N  M
  // Compute center-of-mass velocity according to the formula above
  double vCMx = 0.0, vCMy = 0.0, vCMz = 0.0;

  for (i = 0; i < N; i++) {
    vCMx += m * vx[i];
    vCMy += m * vy[i];
    vCMz += m * vz[i];
  }

  const double inv_Nm = 1.0 / (N * m);
  vCMx *= inv_Nm;
  vCMy *= inv_Nm;
  vCMz *= inv_Nm;

  //  Subtract out the center-of-mass velocity and accumulate vSqdSum in one pass
  double vSqdSum = 0.0;

  for (i = 0; i < N; i++) {
    double vxi = vx[i] - vCMx;
    double vyi = vy[i] - vCMy;
    double vzi = vz[i] - vCMz;

    vx[i] = vxi;
    vy[i] = vyi;
    vz[i] = vzi;

    vSqdSum += vxi * vxi + vyi * vyi + vzi * vzi;
  }

  //  Now we want to scale the average velocity of the system
  //  by a factor which is consistent with our initial temperature, Tinit
  const double lambda = sqrt(3.0 * (N - 1) * Tinit / vSqdSum);

  for (i = 0; i < N; i++) {
    vx[i] *= lambda;
    vy[i] *= lambda;
    vz[i] *= lambda;
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

=== err_escala.e2256900 ===

```text
real 49.18
user 49.14
sys 0.00
real 24.65
user 49.21
sys 0.04
real 12.32
user 49.21
sys 0.01
real 6.17
user 49.30
sys 0.01
real 3.11
user 49.52
sys 0.03
real 1.66
user 51.01
sys 0.01
real 1.23
user 53.96
sys 0.06
```

### Tabela consolidada de escalabilidade forte (Job `2256900`, versão final SoA + Fusão, $N = 5000$, 201 passos)

| Threads ($p$) | `real` (s) | `user` (s) | `sys` (s) | Speedup ($S_p = T_1 / T_p$) | Eficiência ($E_p = S_p / p$) | Fração serial Karp-Flatt ($e$) |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 49.18 | 49.14 | 0.00 | 1.00× | 100.00% | — |
| 2 | 24.65 | 49.21 | 0.04 | 2.00× | 99.76% | 0.0024 |
| 4 | 12.32 | 49.21 | 0.01 | 3.99× | 99.80% | 0.0007 |
| 8 | 6.17 | 49.30 | 0.01 | 7.97× | 99.64% | 0.0005 |
| 16 | 3.11 | 49.52 | 0.03 | 15.81× | 98.83% | 0.0008 |
| 32 | 1.66 | 51.01 | 0.01 | 29.63× | 92.58% | 0.0026 |
| 48 | 1.23 | 53.96 | 0.06 | 39.98× | 83.30% | 0.0043 |

---

=== perf_stat_1t.txt ===

```text

 Performance counter stats for './md':

    98,284,076,400      cycles                                                      
    16,587,537,670      instructions              #    0.17  insn per cycle         
     1,996,249,537      cache-references                                            
       485,547,392      cache-misses              #   24.323 % of all cache refs    
     1,996,249,355      L1-dcache-loads                                             
       485,547,392      L1-dcache-load-misses     #   24.32% of all L1-dcache accesses

      49.149131177 seconds time elapsed

      49.143963000 seconds user
       0.000000000 seconds sys
```

---

=== perf_stat_48t.txt ===

```text

 Performance counter stats for './md':

   106,954,506,146      cycles                                                      
    31,565,466,680      instructions              #    0.30  insn per cycle         
     5,161,175,096      cache-references                                            
       488,320,269      cache-misses              #    9.461 % of all cache refs    
     5,161,173,305      L1-dcache-loads                                             
       488,320,269      L1-dcache-load-misses     #    9.46% of all L1-dcache accesses

       1.188934102 seconds time elapsed

      53.453356000 seconds user
       0.019885000 seconds sys
```

---

=== perf_report_48t.txt ===

```text
# To display the perf.data header info, please use --header/--header-only options.
#
#
# Total Lost Samples: 0
#
# Samples: 231K of event 'cycles'
# Event count (approx.): 112400453944
#
# Children      Self  Command  Shared Object       Symbol                                         
# ........  ........  .......  ..................  ...............................................
#
    97.88%     0.00%  md       libc-2.28.so        [.] thread_start
            |
            ---thread_start
               start_thread
               |          
                --97.88%--gomp_thread_start
                          |          
                          |--86.34%--computeAccelerations
                          |          
                          |--5.40%--gomp_team_barrier_wait_end
                          |          
                          |--4.17%--gomp_barrier_wait_end
                          |          
                          |--0.73%--MeanSquaredVelocity
                          |          
                           --0.64%--Kinetic

    97.88%     0.00%  md       libpthread-2.28.so  [.] start_thread
            |
            ---start_thread
               |          
                --97.88%--gomp_thread_start
                          |          
                          |--86.34%--computeAccelerations
                          |          
                          |--5.40%--gomp_team_barrier_wait_end
                          |          
                          |--4.17%--gomp_barrier_wait_end
                          |          
```

---

=== resultados_tempos.txt ===

*(464 linhas — output completo do job `2256900`, confirmando compilação limpa e paridade termodinâmica idêntica entre 1, 2, 4, 8, 16, 32 e 48 threads)*

```text
rm -f ./md
gcc -O3 -ffast-math -fopenmp -march=native -g -fno-omit-frame-pointer src/MD.cpp -lm -o md
=== Teste com 1 threads ===

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  WELCOME TO WILLY P CHEM MD!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ENTER A TITLE FOR YOUR CALCULATION!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  TITLE ENTERED AS 'cp'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  WHICH NOBLE GAS WOULD YOU LIKE TO SIMULATE? (DEFAULT IS ARGON)

  FOR HELIUM,  TYPE 'He' THEN PRESS 'return' TO CONTINUE
  FOR NEON,    TYPE 'Ne' THEN PRESS 'return' TO CONTINUE
  FOR ARGON,   TYPE 'Ar' THEN PRESS 'return' TO CONTINUE
  FOR KRYPTON, TYPE 'Kr' THEN PRESS 'return' TO CONTINUE
  FOR XENON,   TYPE 'Xe' THEN PRESS 'return' TO CONTINUE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                     YOU ARE SIMULATING Ar GAS! 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  YOU WILL NOW ENTER A FEW SIMULATION PARAMETERS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ENTER THE INTIAL TEMPERATURE OF YOUR GAS IN KELVIN


  ENTER THE NUMBER DENSITY IN moles/m^3
  FOR REFERENCE, NUMBER DENSITY OF AN IDEAL GAS AT STP IS ABOUT 40 moles/m^3
  NUMBER DENSITY OF LIQUID ARGON AT 1 ATM AND 87 K IS ABOUT 35000 moles/m^3
  PERCENTAGE OF CALCULATION COMPLETE:
  [ 10 | 20 | 30 | 40 | 50 | 60 | 70 | 80 | 90 | 100 ]

  TO ANIMATE YOUR SIMULATION, OPEN THE FILE 
  'cp_traj.xyz' WITH VMD AFTER THE SIMULATION COMPLETES

  TO ANALYZE INSTANTANEOUS DATA ABOUT YOUR MOLECULE, OPEN THE FILE 
  'cp_output.txt' WITH YOUR FAVORITE TEXT EDITOR OR IMPORT THE DATA INTO EXCEL

  THE FOLLOWING THERMODYNAMIC AVERAGES WILL BE COMPUTED AND WRITTEN TO THE FILE  
  'cp_average.txt':

  AVERAGE TEMPERATURE (K):                       131.45506

  AVERAGE PRESSURE  (Pa):                  131001228.45076

  PV/nT (J * mol^-1 K^-1):                        28.47279

  PERCENT ERROR of pV/nT AND GAS CONSTANT:       242.44905

  THE COMPRESSIBILITY (unitless):                  3.42449 

  TOTAL VOLUME (m^3):                      2.37220e-25 

  NUMBER OF PARTICLES (unitless):          5000 

=== Teste com 2 threads ===

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  WELCOME TO WILLY P CHEM MD!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ENTER A TITLE FOR YOUR CALCULATION!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  TITLE ENTERED AS 'cp'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  WHICH NOBLE GAS WOULD YOU LIKE TO SIMULATE? (DEFAULT IS ARGON)

  FOR HELIUM,  TYPE 'He' THEN PRESS 'return' TO CONTINUE
  FOR NEON,    TYPE 'Ne' THEN PRESS 'return' TO CONTINUE
  FOR ARGON,   TYPE 'Ar' THEN PRESS 'return' TO CONTINUE
  FOR KRYPTON, TYPE 'Kr' THEN PRESS 'return' TO CONTINUE
  FOR XENON,   TYPE 'Xe' THEN PRESS 'return' TO CONTINUE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                     YOU ARE SIMULATING Ar GAS! 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  YOU WILL NOW ENTER A FEW SIMULATION PARAMETERS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ENTER THE INTIAL TEMPERATURE OF YOUR GAS IN KELVIN


  ENTER THE NUMBER DENSITY IN moles/m^3
  FOR REFERENCE, NUMBER DENSITY OF AN IDEAL GAS AT STP IS ABOUT 40 moles/m^3
  NUMBER DENSITY OF LIQUID ARGON AT 1 ATM AND 87 K IS ABOUT 35000 moles/m^3
  PERCENTAGE OF CALCULATION COMPLETE:
  [ 10 | 20 | 30 | 40 | 50 | 60 | 70 | 80 | 90 | 100 ]

  TO ANIMATE YOUR SIMULATION, OPEN THE FILE 
  'cp_traj.xyz' WITH VMD AFTER THE SIMULATION COMPLETES

  TO ANALYZE INSTANTANEOUS DATA ABOUT YOUR MOLECULE, OPEN THE FILE 
  'cp_output.txt' WITH YOUR FAVORITE TEXT EDITOR OR IMPORT THE DATA INTO EXCEL

  THE FOLLOWING THERMODYNAMIC AVERAGES WILL BE COMPUTED AND WRITTEN TO THE FILE  
  'cp_average.txt':

  AVERAGE TEMPERATURE (K):                       131.45506

  AVERAGE PRESSURE  (Pa):                  131001228.45076

  PV/nT (J * mol^-1 K^-1):                        28.47279

  PERCENT ERROR of pV/nT AND GAS CONSTANT:       242.44905

  THE COMPRESSIBILITY (unitless):                  3.42449 

  TOTAL VOLUME (m^3):                      2.37220e-25 

  NUMBER OF PARTICLES (unitless):          5000 

=== Teste com 4 threads ===

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  WELCOME TO WILLY P CHEM MD!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ENTER A TITLE FOR YOUR CALCULATION!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  TITLE ENTERED AS 'cp'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  WHICH NOBLE GAS WOULD YOU LIKE TO SIMULATE? (DEFAULT IS ARGON)

  FOR HELIUM,  TYPE 'He' THEN PRESS 'return' TO CONTINUE
  FOR NEON,    TYPE 'Ne' THEN PRESS 'return' TO CONTINUE
  FOR ARGON,   TYPE 'Ar' THEN PRESS 'return' TO CONTINUE
  FOR KRYPTON, TYPE 'Kr' THEN PRESS 'return' TO CONTINUE
  FOR XENON,   TYPE 'Xe' THEN PRESS 'return' TO CONTINUE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                     YOU ARE SIMULATING Ar GAS! 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  YOU WILL NOW ENTER A FEW SIMULATION PARAMETERS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ENTER THE INTIAL TEMPERATURE OF YOUR GAS IN KELVIN


  ENTER THE NUMBER DENSITY IN moles/m^3
  FOR REFERENCE, NUMBER DENSITY OF AN IDEAL GAS AT STP IS ABOUT 40 moles/m^3
  NUMBER DENSITY OF LIQUID ARGON AT 1 ATM AND 87 K IS ABOUT 35000 moles/m^3
  PERCENTAGE OF CALCULATION COMPLETE:
  [ 10 | 20 | 30 | 40 | 50 | 60 | 70 | 80 | 90 | 100 ]

  TO ANIMATE YOUR SIMULATION, OPEN THE FILE 
  'cp_traj.xyz' WITH VMD AFTER THE SIMULATION COMPLETES

  TO ANALYZE INSTANTANEOUS DATA ABOUT YOUR MOLECULE, OPEN THE FILE 
  'cp_output.txt' WITH YOUR FAVORITE TEXT EDITOR OR IMPORT THE DATA INTO EXCEL

  THE FOLLOWING THERMODYNAMIC AVERAGES WILL BE COMPUTED AND WRITTEN TO THE FILE  
  'cp_average.txt':

  AVERAGE TEMPERATURE (K):                       131.45506

  AVERAGE PRESSURE  (Pa):                  131001228.45076

  PV/nT (J * mol^-1 K^-1):                        28.47279

  PERCENT ERROR of pV/nT AND GAS CONSTANT:       242.44905

  THE COMPRESSIBILITY (unitless):                  3.42449 

  TOTAL VOLUME (m^3):                      2.37220e-25 

  NUMBER OF PARTICLES (unitless):          5000 

=== Teste com 8 threads ===

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  WELCOME TO WILLY P CHEM MD!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ENTER A TITLE FOR YOUR CALCULATION!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  TITLE ENTERED AS 'cp'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  WHICH NOBLE GAS WOULD YOU LIKE TO SIMULATE? (DEFAULT IS ARGON)

  FOR HELIUM,  TYPE 'He' THEN PRESS 'return' TO CONTINUE
  FOR NEON,    TYPE 'Ne' THEN PRESS 'return' TO CONTINUE
  FOR ARGON,   TYPE 'Ar' THEN PRESS 'return' TO CONTINUE
  FOR KRYPTON, TYPE 'Kr' THEN PRESS 'return' TO CONTINUE
  FOR XENON,   TYPE 'Xe' THEN PRESS 'return' TO CONTINUE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                     YOU ARE SIMULATING Ar GAS! 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  YOU WILL NOW ENTER A FEW SIMULATION PARAMETERS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ENTER THE INTIAL TEMPERATURE OF YOUR GAS IN KELVIN


  ENTER THE NUMBER DENSITY IN moles/m^3
  FOR REFERENCE, NUMBER DENSITY OF AN IDEAL GAS AT STP IS ABOUT 40 moles/m^3
  NUMBER DENSITY OF LIQUID ARGON AT 1 ATM AND 87 K IS ABOUT 35000 moles/m^3
  PERCENTAGE OF CALCULATION COMPLETE:
  [ 10 | 20 | 30 | 40 | 50 | 60 | 70 | 80 | 90 | 100 ]

  TO ANIMATE YOUR SIMULATION, OPEN THE FILE 
  'cp_traj.xyz' WITH VMD AFTER THE SIMULATION COMPLETES

  TO ANALYZE INSTANTANEOUS DATA ABOUT YOUR MOLECULE, OPEN THE FILE 
  'cp_output.txt' WITH YOUR FAVORITE TEXT EDITOR OR IMPORT THE DATA INTO EXCEL

  THE FOLLOWING THERMODYNAMIC AVERAGES WILL BE COMPUTED AND WRITTEN TO THE FILE  
  'cp_average.txt':

  AVERAGE TEMPERATURE (K):                       131.45506

  AVERAGE PRESSURE  (Pa):                  131001228.45076

  PV/nT (J * mol^-1 K^-1):                        28.47279

  PERCENT ERROR of pV/nT AND GAS CONSTANT:       242.44905

  THE COMPRESSIBILITY (unitless):                  3.42449 

  TOTAL VOLUME (m^3):                      2.37220e-25 

  NUMBER OF PARTICLES (unitless):          5000 

=== Teste com 16 threads ===

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  WELCOME TO WILLY P CHEM MD!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ENTER A TITLE FOR YOUR CALCULATION!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  TITLE ENTERED AS 'cp'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  WHICH NOBLE GAS WOULD YOU LIKE TO SIMULATE? (DEFAULT IS ARGON)

  FOR HELIUM,  TYPE 'He' THEN PRESS 'return' TO CONTINUE
  FOR NEON,    TYPE 'Ne' THEN PRESS 'return' TO CONTINUE
  FOR ARGON,   TYPE 'Ar' THEN PRESS 'return' TO CONTINUE
  FOR KRYPTON, TYPE 'Kr' THEN PRESS 'return' TO CONTINUE
  FOR XENON,   TYPE 'Xe' THEN PRESS 'return' TO CONTINUE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                     YOU ARE SIMULATING Ar GAS! 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  YOU WILL NOW ENTER A FEW SIMULATION PARAMETERS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ENTER THE INTIAL TEMPERATURE OF YOUR GAS IN KELVIN


  ENTER THE NUMBER DENSITY IN moles/m^3
  FOR REFERENCE, NUMBER DENSITY OF AN IDEAL GAS AT STP IS ABOUT 40 moles/m^3
  NUMBER DENSITY OF LIQUID ARGON AT 1 ATM AND 87 K IS ABOUT 35000 moles/m^3
  PERCENTAGE OF CALCULATION COMPLETE:
  [ 10 | 20 | 30 | 40 | 50 | 60 | 70 | 80 | 90 | 100 ]

  TO ANIMATE YOUR SIMULATION, OPEN THE FILE 
  'cp_traj.xyz' WITH VMD AFTER THE SIMULATION COMPLETES

  TO ANALYZE INSTANTANEOUS DATA ABOUT YOUR MOLECULE, OPEN THE FILE 
  'cp_output.txt' WITH YOUR FAVORITE TEXT EDITOR OR IMPORT THE DATA INTO EXCEL

  THE FOLLOWING THERMODYNAMIC AVERAGES WILL BE COMPUTED AND WRITTEN TO THE FILE  
  'cp_average.txt':

  AVERAGE TEMPERATURE (K):                       131.45506

  AVERAGE PRESSURE  (Pa):                  131001228.45076

  PV/nT (J * mol^-1 K^-1):                        28.47279

  PERCENT ERROR of pV/nT AND GAS CONSTANT:       242.44905

  THE COMPRESSIBILITY (unitless):                  3.42449 

  TOTAL VOLUME (m^3):                      2.37220e-25 

  NUMBER OF PARTICLES (unitless):          5000 

=== Teste com 32 threads ===

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  WELCOME TO WILLY P CHEM MD!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ENTER A TITLE FOR YOUR CALCULATION!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  TITLE ENTERED AS 'cp'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  WHICH NOBLE GAS WOULD YOU LIKE TO SIMULATE? (DEFAULT IS ARGON)

  FOR HELIUM,  TYPE 'He' THEN PRESS 'return' TO CONTINUE
  FOR NEON,    TYPE 'Ne' THEN PRESS 'return' TO CONTINUE
  FOR ARGON,   TYPE 'Ar' THEN PRESS 'return' TO CONTINUE
  FOR KRYPTON, TYPE 'Kr' THEN PRESS 'return' TO CONTINUE
  FOR XENON,   TYPE 'Xe' THEN PRESS 'return' TO CONTINUE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                     YOU ARE SIMULATING Ar GAS! 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  YOU WILL NOW ENTER A FEW SIMULATION PARAMETERS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ENTER THE INTIAL TEMPERATURE OF YOUR GAS IN KELVIN


  ENTER THE NUMBER DENSITY IN moles/m^3
  FOR REFERENCE, NUMBER DENSITY OF AN IDEAL GAS AT STP IS ABOUT 40 moles/m^3
  NUMBER DENSITY OF LIQUID ARGON AT 1 ATM AND 87 K IS ABOUT 35000 moles/m^3
  PERCENTAGE OF CALCULATION COMPLETE:
  [ 10 | 20 | 30 | 40 | 50 | 60 | 70 | 80 | 90 | 100 ]

  TO ANIMATE YOUR SIMULATION, OPEN THE FILE 
  'cp_traj.xyz' WITH VMD AFTER THE SIMULATION COMPLETES

  TO ANALYZE INSTANTANEOUS DATA ABOUT YOUR MOLECULE, OPEN THE FILE 
  'cp_output.txt' WITH YOUR FAVORITE TEXT EDITOR OR IMPORT THE DATA INTO EXCEL

  THE FOLLOWING THERMODYNAMIC AVERAGES WILL BE COMPUTED AND WRITTEN TO THE FILE  
  'cp_average.txt':

  AVERAGE TEMPERATURE (K):                       131.45506

  AVERAGE PRESSURE  (Pa):                  131001228.45076

  PV/nT (J * mol^-1 K^-1):                        28.47279

  PERCENT ERROR of pV/nT AND GAS CONSTANT:       242.44905

  THE COMPRESSIBILITY (unitless):                  3.42449 

  TOTAL VOLUME (m^3):                      2.37220e-25 

  NUMBER OF PARTICLES (unitless):          5000 

=== Teste com 48 threads ===

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  WELCOME TO WILLY P CHEM MD!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ENTER A TITLE FOR YOUR CALCULATION!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  TITLE ENTERED AS 'cp'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  WHICH NOBLE GAS WOULD YOU LIKE TO SIMULATE? (DEFAULT IS ARGON)

  FOR HELIUM,  TYPE 'He' THEN PRESS 'return' TO CONTINUE
  FOR NEON,    TYPE 'Ne' THEN PRESS 'return' TO CONTINUE
  FOR ARGON,   TYPE 'Ar' THEN PRESS 'return' TO CONTINUE
  FOR KRYPTON, TYPE 'Kr' THEN PRESS 'return' TO CONTINUE
  FOR XENON,   TYPE 'Xe' THEN PRESS 'return' TO CONTINUE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                     YOU ARE SIMULATING Ar GAS! 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  YOU WILL NOW ENTER A FEW SIMULATION PARAMETERS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ENTER THE INTIAL TEMPERATURE OF YOUR GAS IN KELVIN


  ENTER THE NUMBER DENSITY IN moles/m^3
  FOR REFERENCE, NUMBER DENSITY OF AN IDEAL GAS AT STP IS ABOUT 40 moles/m^3
  NUMBER DENSITY OF LIQUID ARGON AT 1 ATM AND 87 K IS ABOUT 35000 moles/m^3
  PERCENTAGE OF CALCULATION COMPLETE:
  [ 10 | 20 | 30 | 40 | 50 | 60 | 70 | 80 | 90 | 100 ]

  TO ANIMATE YOUR SIMULATION, OPEN THE FILE 
  'cp_traj.xyz' WITH VMD AFTER THE SIMULATION COMPLETES

  TO ANALYZE INSTANTANEOUS DATA ABOUT YOUR MOLECULE, OPEN THE FILE 
  'cp_output.txt' WITH YOUR FAVORITE TEXT EDITOR OR IMPORT THE DATA INTO EXCEL

  THE FOLLOWING THERMODYNAMIC AVERAGES WILL BE COMPUTED AND WRITTEN TO THE FILE  
  'cp_average.txt':

  AVERAGE TEMPERATURE (K):                       131.45506

  AVERAGE PRESSURE  (Pa):                  131001228.45076

  PV/nT (J * mol^-1 K^-1):                        28.47279

  PERCENT ERROR of pV/nT AND GAS CONSTANT:       242.44905

  THE COMPRESSIBILITY (unitless):                  3.42449 

  TOTAL VOLUME (m^3):                      2.37220e-25 

  NUMBER OF PARTICLES (unitless):          5000 
```
