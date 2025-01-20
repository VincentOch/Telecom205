# Spécifications amplificateurs

| Modèle d'ampli | Gain | OIP3 | NF | Pc@1dB | Tension | Consommation |
| --- | --- | --- | --- | --- | --- | --- |
| ZX60-V63+ | 20 dB | 30 dBm | 3.7 | 17.8 dBm | 5 V | 69 mA |
| ZX60-V62+ | 15.4 dB | 33 dBm | 5.1 | 18.5 dBm | 5 V | 82 mA |
| ZHL-42    | 33 dB | 38 dBm | 7.55 | 28.32 dBm | 15 V | 880 mA |
| RFLUPA05M06G | 33 dB | 40.5 dBm | 3 | 30 dBm | 12 V | 280 mA |
| ADL5606   | 24.3 dB | 45.5 dBm | 4.7 | 30.8 dBm | 5 V | 362 mA |

> Pour déterminer le IPP3 : $$IPP3 = P_{in} + \frac{IM_3}{2}$$
> Pour déterminer le OIP3 : $$OIP3 = P_{out} + \frac{IM_3}{2}$$
> Pour déterminer le OIP3 à partir de l'IPP3 (en prenant un point en régime linaire, loin de la compression) : $$OIP3 = OIP_3 - Gain$$

On peut changer le Vref du DAC ou mettre un gain négatif avant l'ampli
