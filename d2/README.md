# Environnement de développement des codes pour TELECOM20x

## TELECOM201b

Dans TELECOM201b, on doit utiliser des codes mettant en oeuvre du filtrage analogique. 

Les scripts principaux utilisé dans le TP sont : 

| File script                         | Involves analog filter |
|-------------------------------------|------------------------|
| `DAC_TB.m`                          |                        |
| `ButterworthFilterSpecifications.m` |        &check;         |
| `TX_BasebandChain.m`                |        &check;         |
| `TX_Chain_woPA.m`                   |        &check;         |
| `TX_Chain.m`                        |        &check;         |
| `LNA_TB.m`                          |                        |
| `RX_TB.m`                           |        &check;         |
| `Play_OriginalAudio.m`              |                        |
| `TXRX.m`                            |        &check;         |

Il faut noter que les fonctions (dans `subblocks/`) suivantes utilisent du filtrage analogique:

| File script                         | Involves analog filter |
|-------------------------------------|------------------------|
| `subblocks/RX.m`                    |        &check;         |
| `subblocks/TX.m`                    |        &check;         |
| `subblocks/basebandAnalogFilt.m`    |        &check;         |








## TELECOM205

Dans TELECOM205, on choisit d'implémeter le filtrage par des filtres FIR. 
On fait une version `_proj` pour les scripts et fonctions utiles. 

| TELECOM201b                         | TELECOM205                            |
|-------------------------------------|---------------------------------------|
| `DAC_TB.m`                          |                                       |
| `ButterworthFilterSpecifications.m` |  N/A                                  |
| `TX_BasebandChain.m`                |  `TX_BasebandChain_proj.m`            |
| `TX_Chain_woPA.m`                   |  `TX_Chain_woPA_proj.m`               |
| `TX_Chain.m`                        |  `TX_Chain_proj.m`                    |
| `LNA_TB.m`                          |  N/A                                  |
| `RX_TB.m`                           |  `RX_TB_proj.m`                       |
| `Play_OriginalAudio.m`              |  N/A                                  |
| `TXRX.m`                            |  `TXRX_proj.m`                        |
|  N/A                                |  `completeTxRx_proj.m`                |
|  N/A                                |  `TX_TB_proj.m`                       |
| `subblocks/RX.m`                    |  `subblocks/RX_proj.m`                |
| `subblocks/TX.m`                    |  `subblocks/TX_proj.m`                |
| `subblocks/basebandAnalogFilt.m`    |  `subblocks/basebandAnalogFiltFake.m` |



Il faut noter que les fonctions (dans `subblocks/`) suivantes utilisent du filtrage FIR:

- `subblocks/BBamp.m`



## Makefile

### TELECOM205

Le makefile doit rassembler les fichiers du projet et les placer dans une archive qui va sur: 

```
c2s.telecom-paristech.fr:TELECOM201/documents/projet/
```

### TELECOM201

Il n'y a pas de Makefile pour TELECOM201 ici car, ils sont implémentés dans: 

```
https://gitlab.telecom-paris.fr/telecom201/labs/frontend-rf-telecom201
https://gitlab.telecom-paris.fr/telecom201/ics905/frontend-rf-telecom201
https://gitlab.telecom-paris.fr/telecom201/labs/tp_filteradc
```
