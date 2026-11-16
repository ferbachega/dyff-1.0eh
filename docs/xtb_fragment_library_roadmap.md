# Biblioteca de fragmentos xTB para refinamento de parâmetros (DYFF → OPLS)

*Roteiro de implementação — EasyHybrid 3, 2026-10-06. Nada implementado ainda.*

## 1. Objetivo

Uma ferramenta que:

1. identifica, para cada termo ligado de uma molécula (ligação, ângulo, diedro), o **subsistema mínimo** (fragmento) que reproduz o ambiente químico daquele termo;
2. faz a varredura desse fragmento com **xTB** (diedros: varredura relaxada; ângulos e ligações: Hessiana + varreduras curtas de validação);
3. guarda cada resultado numa **base de dados** indexada por um **código de ambiente químico**, de modo que o mesmo ambiente nunca seja calculado duas vezes;
4. usa os dados para **refinar os parâmetros do DYFF** (depois OPLS) de qualquer molécula que contenha ambientes já conhecidos;
5. torna trivial **adicionar ambientes desconhecidos** (triagem automática a partir de qualquer molécula, SMILES, arquivo ou seleção).

O que já existe e será reaproveitado:

| Peça | Onde está |
|---|---|
| Varredura relaxada/rígida xTB, ajuste linear de séries de cossenos, curva verificada pelo pDynamo | `src/pdynamo/torsion_fit.py` |
| Edição de parâmetros DYFF sem tocar no pDynamo (wrapper das *factories*) | `src/pdynamo/dyff_edits.py` |
| Construção de sistemas DYFF a partir de objetos do Builder; cargas CM5 | `src/gui/windows/builder/dyff_parameters.py`, `src/pdynamo/opls/ligand.py` |
| Janela de ajuste com EasyPlot | `src/gui/windows/builder/torsion_fit_window.py` |
| 663 moléculas curadas (biomoléculas, fármacos, solventes…) | `src/gui/windows/builder/structures/` |
| RDKit 2025.09 (Python), OpenBabel 3.1 (CLI), xtb 6.7.1, SQLite 3.37 | ambiente atual |

---

## 2. Arquitetura geral

```
molécula (Builder / sistema / SMILES / arquivo)
        │
        ▼
[A] enumeração de termos  ──►  [B] código de ambiente (por termo)
                                         │
                                         ▼
                               [C] consulta à base de dados
                         ┌──────────────┴───────────────┐
                    conhecido                      desconhecido
                         │                              │
                         ▼                              ▼
              [F] ajuste MM no fragmento     [D] fragmentação mínima
                  (milissegundos)                       │
                         │                              ▼
                         │                   [E] fila de cálculos xTB
                         │                     (varredura / Hessiana)
                         │                              │
                         │                              ▼
                         │                     grava na base de dados
                         ▼                              │
              [G] aplica ao DYFF  ◄─────────────────────┘
              (edição por termo, Undo, relatório de cobertura)
```

Princípio central: **a base guarda dados de mecânica quântica (QM), não parâmetros**. A varredura xTB de um fragmento é transferível e cara; o ajuste do lado MM é barato e depende do modelo MM corrente (cargas, Lennard-Jones, demais termos). Portanto o ajuste é refeito, ou lido de um cache secundário, sempre que o modelo MM muda. Assim a mesma varredura serve ao DYFF hoje e ao OPLS depois, sem recalcular nada com xTB.

---

## 3. Código do ambiente químico

O código precisa ser **determinístico** (a mesma química dá o mesmo código, independente da numeração dos átomos e do arquivo de origem) e **hierárquico** (para permitir reaproveitar dados de ambientes parecidos quando não há um idêntico).

### 3.1 Dois níveis de identidade

**(a) Identidade exata do fragmento (chave do cálculo):**

```
FRAG = SMILES canônico do fragmento com capeamento (RDKit, com estereo E/Z quando relevante)
TERM = tipo + índices canônicos dos átomos do termo no fragmento
       ex.: "dihedral:4-7-8-11"  (rank canônico do RDKit, orientação normalizada)
CHAVE_EXATA = sha1( FRAG | TERM | carga | multiplicidade )[:16]
```

Dois usuários que importarem a mesma molécula de arquivos diferentes geram a mesma chave. O SMILES vai junto, legível, para inspeção.

**(b) Código de ambiente (chave de transferência), em raios crescentes:**

Para cada termo, com os átomos do termo como centro, gera-se um código no estilo Morgan/ECFP restrito à vizinhança:

| Nível | Conteúdo | Uso |
|---|---|---|
| `r0` | tipo do termo + tipos DYFF dos átomos + ordens das ligações do termo (a chave DYFF atual) | fallback mínimo (= hoje) |
| `r1` | `r0` + elemento/hibridização/aromaticidade/carga/nº de H dos vizinhos diretos dos átomos do termo | transferência boa para ligações e ângulos |
| `r2` | `r1` + mais uma camada (vizinhos dos vizinhos), anéis marcados por tamanho | transferência para diedros (efeitos 1-5, conjugação) |
| `exato` | CHAVE_EXATA | reaproveitamento sem perda |

Formato legível + hash, por exemplo:

```
dihedral|r2|C.ar:C.ar-C.3:O.3|9f3a1c0e7b2d4a61
```

Consulta: **exato → r2 → r1 → r0 (parâmetro original do DYFF)**. O nível usado entra no relatório de cobertura, com uma cor por nível.

Implementação: RDKit (`Chem.MolToSmiles` com mapas de átomos para os índices canônicos; `rdFingerprintGenerator`/`GetMorganGenerator` com `fromAtoms` = átomos do termo para os hashes de ambiente). Os invariantes de átomo são definidos por nós (elemento, grau, H, carga, aromaticidade, anel), não os padrões do RDKit, para ficarem estáveis entre versões.

### 3.2 Por que não usar só a chave DYFF

A chave DYFF (`O:Res:C:Res:C:Res:C:Res:S`) é genérica demais. Um diedro ajustado na flavona seria aplicado a qualquer arila ligada a oxigênio. Os níveis `r1`/`r2` resolvem isso, e a aplicação passa a ser **por termo** (seção 7), não por chave.

---

## 4. Fragmentação: subsistema mínimo

Referência metodológica: *Stern, C. D. et al., "Capturing non-local through-bond effects in molecular mechanics force fields I: Fragmenting molecules for quantum chemical torsion scans", J. Chem. Theory Comput. (2020)*, o "fragmenter" do Open Force Field. Usaremos a mesma ideia de validar o fragmento pela **ordem de ligação de Wiberg (WBO)**, que o xtb já fornece de graça (arquivo `wbo`).

### 4.1 Regras de construção (para um termo com átomos T)

1. **Núcleo:** átomos de T + vizinhos diretos (para diedro: os átomos 1-4 e seus vizinhos).
2. **Não cortar:**
   - anéis que contenham qualquer átomo do núcleo (o anel entra inteiro);
   - sistemas conjugados ligados ao núcleo (ligações duplas/aromáticas contíguas), até um limite de tamanho;
   - grupos funcionais "atômicos": amida, carboxila/carboxilato, éster, nitro, sulfonila, fosfato, guanidínio, imidazol etc. (lista por SMARTS, extensível pelo usuário).
3. **Cortar só em ligações simples, preferencialmente C(sp³)–C(sp³)**, a pelo menos 2 ligações dos átomos de T.
4. **Capeamento:** H no lugar do átomo removido, na direção original, com comprimento X–H padrão. Opção de capear com CH₃ quando o átomo cortado é heteroátomo ou carbono ligado a grupo retirador, para preservar o efeito indutivo.
5. **Carga e multiplicidade:** soma das cargas formais dos átomos mantidos; recusa radicais acidentais.

### 4.2 Validação pela WBO (só para diedros)

- Calcula a WBO da ligação central na **molécula inteira** (um *single point* xtb, segundos) e no fragmento.
- Se |ΔWBO| > 0,03, **cresce o fragmento** uma camada no lado que mais contribui e repete. É o critério do OpenFF.
- Limite de tamanho, por exemplo 40 átomos pesados. Se ultrapassar, usa o maior fragmento válido e marca o termo como "fragmento truncado".

### 4.3 Deduplicação

Antes de enfileirar qualquer cálculo, todos os fragmentos de todas as moléculas da triagem passam pela CHAVE_EXATA. Uma proteína com 50 resíduos de serina gera **um** fragmento por tipo de termo, não 50.

### 4.4 Custo estimado (xtb GFN2, 1 CPU, fragmentos de 15–30 átomos)

| Cálculo | Tempo típico |
|---|---|
| WBO / *single point* | < 1 s |
| Otimização + Hessiana (todas as ligações e ângulos do fragmento) | 2–10 s |
| Varredura relaxada de diedro (24 pontos, passo 15°) + *single points* | 5–30 s |
| 1 000 fragmentos com diedro | ~2–8 h em 1 núcleo; minutos em paralelo |

---

## 5. O que é calculado por tipo de termo

| Termo | Cálculo principal | Validação opcional | Parâmetros DYFF obtidos |
|---|---|---|---|
| **Diedro** | varredura relaxada 360° (`torsion_fit.run_scan`), sentido duplo opcional | varredura rígida | série Σ cₙ cos(nφ), já suportada |
| **Ângulo** | Hessiana no mínimo → **método de Seminario modificado** (*Allen, Dodda, Cole, JCTC 2018*) | varredura de ±25° com restrição (8–10 pontos) | θ₀ e k → forma cosseno UFF/DYFF: c₂ = 1/(4 sin²θ₀), c₁ = −4c₂cos θ₀, c₀ = c₂(2cos²θ₀+1) |
| **Ligação** | Hessiana → Seminario (k), geometria otimizada (b₀) | varredura de ±0,12 Å (7 pontos) | b₀, k (DYFF: k(b−b₀)², sem ½) |
| **Impróprio** | Hessiana (projeção fora do plano) | varredura curta | c₀ do termo fora do plano |

Por que a Hessiana para ligações e ângulos: **uma** Hessiana por fragmento dá as constantes de *todas* as ligações e ângulos dele, contra dezenas de varreduras 1D. As varreduras curtas ficam como verificação, e como modo explícito quando o usuário quiser a curva no EasyPlot, como já fazemos com diedros.

---

## 6. Base de dados

### 6.1 Formato: SQLite

**Recomendação: um arquivo SQLite** (`sqlite3` da biblioteca padrão do Python):

- é um único arquivo, sem servidor, transacional (não corrompe se o EasyHybrid fechar no meio de uma fila) e fácil de copiar e compartilhar;
- aceita consultas por qualquer nível do código de ambiente (índices SQL);
- os dados volumosos (geometrias da varredura) vão como BLOB comprimido (`numpy.savez_compressed`).

Tamanho estimado de uma varredura de diedro: 24 pontos × 30 átomos × 3 × float32 ≈ 9 kB comprimidos. Cem mil varreduras dão cerca de 1 GB, aceitável. As geometrias **precisam ser guardadas**: é o que permite refazer o ajuste MM quando as cargas ou o campo de força mudam, sem xTB.

### 6.2 Duas camadas

```
src/data/fragment_library/seed.sqlite        ← distribuída com o EasyHybrid (somente leitura)
~/.config/EasyHybrid/fragment_library.sqlite ← do usuário (cálculos novos, ambientes próprios)
```

As consultas olham primeiro a base do usuário, depois a semente. A exportação/importação de pacotes (`.ehfrag.zip`: SQLite parcial + manifesto JSON) permite compartilhar entre colegas e grupos. A mescla é feita pela CHAVE_EXATA, guardando o registro de método mais alto, ou o mais recente para o mesmo método.

### 6.3 Esquema (proposta)

```sql
fragments (
  id INTEGER PRIMARY KEY,
  exact_smiles TEXT NOT NULL,         -- canônico, com capeamento
  inchikey TEXT,                      -- conferência cruzada
  charge INTEGER, multiplicity INTEGER,
  n_atoms INTEGER, n_heavy INTEGER,
  xyz_opt BLOB,                       -- geometria otimizada xTB
  cm5 BLOB,                           -- cargas CM5 do fragmento (para o ajuste MM)
  source TEXT,                        -- "seed:structures/26_cofactors/NAD", "user:smiles", "ccd:ATP"...
  created TEXT, UNIQUE(exact_smiles, charge, multiplicity)
);

terms (
  id INTEGER PRIMARY KEY,
  fragment_id INTEGER REFERENCES fragments,
  kind TEXT,                          -- bond | angle | dihedral | improper
  atoms TEXT,                         -- índices canônicos "4-7-8-11"
  exact_key TEXT UNIQUE,              -- CHAVE_EXATA
  env_r2 TEXT, env_r1 TEXT, env_r0 TEXT,
  wbo_parent REAL, wbo_fragment REAL, truncated INTEGER
);
CREATE INDEX terms_r2 ON terms(env_r2);  -- e r1, r0

calculations (
  id INTEGER PRIMARY KEY,
  term_id INTEGER REFERENCES terms,     -- ou fragment_id para Hessianas
  kind TEXT,                            -- relaxed_scan | rigid_scan | hessian | short_scan
  method TEXT,                          -- "GFN2-xTB", "GFN2-xTB/ALPB(water)"
  protocol TEXT,                        -- JSON: passo, faixa, sentidos, restrição
  program_version TEXT,                 -- "xtb 6.7.1"
  coordinate BLOB, energies BLOB, geometries BLOB,
  hessian BLOB,
  status TEXT, wall_time REAL, created TEXT
);

fits (                                  -- cache; pode ser apagado e refeito
  id INTEGER PRIMARY KEY,
  calculation_id INTEGER REFERENCES calculations,
  force_field TEXT, ff_version TEXT,    -- "DYFF", "dyff-1.0" / "OPLS", "oplsaa08"
  charge_model TEXT,                    -- "CM5x1.00"
  parameters TEXT,                      -- JSON nativo do campo de força
  rmse REAL, weighted_rmse REAL, created TEXT
);

user_groups (                           -- grupos químicos definidos pelo usuário
  id INTEGER PRIMARY KEY, name TEXT, smarts TEXT, keep_whole INTEGER, example_smiles TEXT
);

meta (key TEXT PRIMARY KEY, value TEXT) -- versão do esquema, versão do código de ambiente
```

O `meta` guarda a **versão do algoritmo de código de ambiente**. Se as regras de hash mudarem, os códigos `r1`/`r2` podem ser recalculados a partir de `exact_smiles` sem refazer nenhum cálculo.

---

## 7. Aplicação ao DYFF

### 7.1 Edição por termo, não por chave

As *factories* do DYFF recebem `(types, indices, connectivity, key)`, ou seja, **os índices dos átomos de cada termo**. O `MMModelDYFFEdited` atual compara só a chave. A extensão é aceitar também edições indexadas pelos átomos:

```python
parameterEdits = {
  "dihedral": { "O:Res:C:Res:C:Res:C:Res:S": [...],        # por chave (hoje)
                (3, 5, 10, 14): [...] },                    # por termo (novo)
  ...
}
```

A ordem de prioridade é: edição por termo → edição por chave → regra do DYFF. Continua picklável e continua sem mudar nada no pDynamo.

### 7.2 Fluxo "Refine from library" (janela DYFF Parameters)

1. enumera os termos da molécula e calcula os códigos;
2. consulta a base (exato → r2 → r1);
3. para cada termo coberto, monta o fragmento a partir da base (geometria otimizada, cargas CM5 do fragmento) e **refaz o ajuste MM com o DYFF atual**, ou lê do cache `fits`;
4. mostra a prévia: termo, nível de correspondência, parâmetro antigo e novo, RMSE do fragmento;
5. "Apply" grava tudo como **uma etapa de Undo** (`dyff_parameter_edits`, como hoje). Passa a valer no Optimize, no Force Field ▸ DYFF e na incorporação ao fechar o Builder.

### 7.3 Redundância de diedros na mesma ligação

Como discutido antes, uma varredura só determina a soma das contribuições dos diedros de uma ligação. Política: **para cada ligação girável, um diedro "portador"** (o de átomos pesados com menor índice canônico) recebe a série ajustada, e os demais diedros daquela ligação são zerados. Assim nada é contado duas vezes e o resultado não depende da ordem do ajuste. É o que o ffTK e o OpenFF BespokeFit fazem, cada um a seu modo.

### 7.4 Cores na janela

Nova cor de fonte "library" (ex.: azul-petróleo) por nível: exato, r2 ou r1. A coluna Quality mostra, por exemplo, `library r2 · 0.4 kJ/mol`. Os termos sem dado ficam com a cor de hoje (força do DYFF) e alimentam a triagem.

---

## 8. Triagem e importação de ambientes desconhecidos

### 8.1 Entradas aceitas

| Entrada | Como |
|---|---|
| Molécula do Builder ou sistema aberto | botão "Library coverage…" na janela DYFF Parameters; menu *Extras ▸ Fragment Library* |
| Seleção (viewing selection) | só os termos dentro da seleção |
| SMILES colado (um por linha) | caixa de texto na janela da biblioteca |
| Arquivos `.mol2 / .sdf / .pdb / .xyz` ou pasta inteira | importação em lote (OpenBabel CLI para converter, RDKit para perceber) |
| Grupo químico novo | nome + SMARTS + molécula exemplo → tabela `user_groups`; passa a ser respeitado pela fragmentação ("não cortar") e vira uma linha de cobertura própria |
| Linha de comando | `python3 -m pdynamo.fragment_library triage moléculas.sdf --run` |

### 8.2 Janela "Fragment Library" (proposta)

```
┌ Fragment Library ─────────────────────────────────────────────────────────┐
│ Source: [Builder molecule ▾]   [Import SMILES/files…]  [Add chemical group…]│
│ Coverage: 61 terms · 44 exact · 9 r2 · 3 r1 · 5 missing (3 fragments)      │
│ ┌────────────┬──────────┬──────────────────────┬───────┬──────────┬──────┐ │
│ │ Term       │ Level    │ Fragment (SMILES)    │ Atoms │ Est. cost│Status│ │
│ │ C9-C4-C2-O1│ missing  │ c1ccc(-c2cc(=O)...   │  27   │  20 s    │queued│ │
│ │ ...                                                                    │ │
│ └────────────────────────────────────────────────────────────────────────┘ │
│ Method [GFN2 ▾] Solvent [none ▾] Parallel jobs [4]  [Run missing] [Cancel] │
│ [3D preview of the selected fragment, cut bonds and caps highlighted]      │
└────────────────────────────────────────────────────────────────────────────┘
```

- **Prévia 3D do fragmento** (reaproveitando o preview da biblioteca de estruturas), com as ligações cortadas e os capeamentos destacados, para o usuário conferir o subsistema.
- **Fila** pelo Process Manager existente: execução em segundo plano, cancelável e retomável. Um cálculo interrompido fica `status = 'pending'` no SQLite e recomeça na próxima sessão.
- Execução paralela de N processos xtb (`OMP_NUM_THREADS` dividido entre os processos).

---

## 9. Fonte para a base semente

Recomendação em três camadas, da mais imediata para a mais ambiciosa.

### 9.1 Camada 1, imediata: nossa própria biblioteca de estruturas
As 663 moléculas de `src/gui/windows/builder/structures/` (aminoácidos, nucleotídeos, açúcares, lipídeos, cofatores, vitaminas, neurotransmissores, fármacos, solventes) já estão curadas e no repositório. Elas cobrem exatamente o que o EasyHybrid simula (biomoléculas e ligantes). Fragmentadas e deduplicadas, devem dar **alguns milhares de fragmentos únicos**. O número exato sai da fase 1, antes de qualquer cálculo xTB.

### 9.2 Camada 2: ligantes do PDB (Chemical Component Dictionary, RCSB)
O CCD tem cerca de 45 mil componentes, com coordenadas ideais e ligações definidas, e já o usamos para os metais. Basta filtrar os componentes que aparecem em mais estruturas do PDB (top 2–5 mil): são os ligantes que os usuários do EasyHybrid realmente vão encontrar. Licença: domínio público (wwPDB).

### 9.3 Camada 3, referência de alto nível: conjuntos do Open Force Field no QCArchive
Os conjuntos "OpenFF Gen2 Torsion Set", "OpenFF Industry Benchmark Torsion Set" e afins contêm **milhares de varreduras de torsão já fragmentadas, em B3LYP-D3BJ/DZVP**. Usos:

- **validar o xTB** nos mesmos fragmentos (seção 10), o que é essencial para uma publicação;
- opcionalmente, guardar essas curvas como registros de método superior (`method = "B3LYP-D3BJ/DZVP"`), usados no lugar do xTB quando existirem.

Licença: os dados do OpenFF/QCArchive são distribuídos sob CC-BY 4.0. **É preciso confirmar por conjunto antes de redistribuir** na base semente; a alternativa é baixá-los sob demanda.

Fontes complementares, se necessário: ChEMBL (fármacos aprovados, para estatística de cobertura) e a lista de grupos funcionais de Ertl (*J. Cheminform.* 2017, implementada no RDKit Contrib), para a tabela inicial de grupos "não cortar".

---

## 10. Validação (necessária para publicar)

1. **Fragmento vs. molécula inteira:** em ~50 moléculas, varredura xTB do diedro na molécula inteira contra a do fragmento (RMSE do perfil, diferença de barreira). Mostra que a fragmentação com WBO preserva o perfil.
2. **xTB vs. DFT:** nas varreduras do OpenFF, RMSE do GFN2-xTB contra B3LYP-D3BJ. Define em quais ambientes o xTB é confiável e marca os outros (ex.: certos heteroaromáticos, conjugação estendida).
3. **DYFF original vs. DYFF refinado:** em moléculas fora da base (conjunto de teste), comparar com o xTB as barreiras de torsão, as geometrias otimizadas (RMSD) e as frequências (Hessiana MM vs. xTB).
4. **Transferência por nível:** erro quando o dado vem de `exato`, `r2` ou `r1`. Justifica a hierarquia e os limiares.

---

## 11. Fases de implementação

| Fase | Conteúdo | Entregável verificável |
|---|---|---|
| **0. Código de ambiente** | módulo `src/pdynamo/fragment_library/environment.py` (RDKit): enumeração de termos, CHAVE_EXATA, hashes r0–r2, ponte objeto do Builder/sistema pDynamo ↔ RDKit | testes: mesma molécula com outra numeração → mesmos códigos; estatística de ambientes únicos nas 663 moléculas |
| **1. Fragmentador** | `fragmenter.py`: regras da seção 4, capeamento, lista SMARTS de grupos "não cortar", validação WBO | testes de fragmentos esperados (flavona, adenosina, ATP, peptídeo); relatório de quantos fragmentos únicos a biblioteca gera |
| **2. Banco de dados** | `store.py`: esquema SQLite, duas camadas, import/export `.ehfrag.zip`, recálculo de códigos | testes de ida e volta; mescla de duas bases |
| **3. Motor de cálculo** | `engine.py`: fila, paralelismo, retomada; varreduras de diedro (reuso de `torsion_fit`), Hessiana + Seminario, varreduras curtas | testes com xtb em 3–4 fragmentos; Seminario conferido contra a Hessiana numérica |
| **4. Ajuste e aplicação ao DYFF** | edição por termo em `dyff_edits.py`; ajuste MM no fragmento com cache `fits`; diedro portador por ligação; botão "Refine from library" | energia com edição por termo = reconstrução real; Undo; flavona refinada vs. xTB |
| **5. Janela e triagem** | janela "Fragment Library" (glade), cobertura, prévia 3D, importação SMILES/arquivos, grupos do usuário, CLI | teste de interface isolado como os anteriores |
| **6. Base semente** | processar a camada 1 (horas de CPU), depois a camada 2; publicar `seed.sqlite` | estatística de cobertura de uma lista de ligantes do PDB |
| **7. Validação** | seção 10 + documento de metodologia (`docs/`) | tabelas/figuras para o artigo |
| **8. OPLS** | adaptador de ajuste para a forma de Fourier (já existe em `opls/torsion_fit.py`) e mapeamento dos termos OPLS por índice | mesmo fluxo na janela OPLS Parameters |

As fases 0–2 não usam xTB nenhuma vez e já respondem à pergunta mais importante antes de gastar CPU: **quantos ambientes únicos existem de fato**.

---

## 12. Riscos e decisões em aberto

| Risco ou decisão | Proposta |
|---|---|
| Precisão do xTB em torsões (heterociclos, conjugação longa) | validação contra o OpenFF (seção 10.2); usar dados DFT quando existirem; marcar ambientes de baixa confiança |
| Código muito específico (pouco reaproveitamento) × muito geral (transferência ruim) | hierarquia exato/r2/r1 com o nível mostrado na interface; limiares decididos pela validação 10.4 |
| Fragmentos carregados em fase gasosa (carboxilatos, aminas protonadas) | opção de ALPB(água) por fragmento carregado; o método fica gravado no registro |
| Mudança de versão do xtb | `program_version` em cada cálculo; nunca misturar automaticamente versões no mesmo ajuste |
| Crescimento da base | geometrias comprimidas; `fits` é um cache descartável; exportar e limpar por fonte |
| Metais e elementos sem DYFF/xTB | fora do escopo inicial (o fragmentador recusa e reporta) |
| Cargas no ajuste MM | CM5 do próprio fragmento (consistente com o protocolo dos ligantes OPLS); gravadas em `fragments.cm5` |
| Estereoquímica e tautômeros | SMILES isomérico na chave exata; tautômeros tratados como ambientes distintos (decisão do usuário na importação) |

---

## 13. Resumo das recomendações

- **Fonte da base:** começar pela biblioteca de estruturas do próprio EasyHybrid (663 moléculas, imediata), depois os ligantes mais frequentes do PDB/CCD, e usar os conjuntos de torsão do Open Force Field (QCArchive, DFT) para validar o xTB.
- **Importação de desconhecidos:** a partir de qualquer molécula aberta ("Library coverage…"), SMILES colado, arquivos ou pasta, seleção, ou um grupo novo definido por SMARTS. A fragmentação automática, validada pela WBO, define o subsistema mínimo e deduplica antes de calcular.
- **Armazenamento:** SQLite em duas camadas (semente distribuída + base do usuário), com dados QM (geometrias e energias comprimidas) como registro primário e parâmetros ajustados como cache recalculável; pacotes `.ehfrag.zip` para compartilhar.
