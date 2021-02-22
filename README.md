# Single_cell_csf
single_cell_csf_last.R - код для кластеризации клеток образцов scRNA-seq CSF и PBMC доноров с MS и контрольных доноров.
Данные взяты из статьи Integrated single cell analysis of blood and cerebrospinal fluid leukocytes in multiple sclerosis (https://doi.org/10.1038/s41467-019-14118-w)
На вход подавались по 6 повторностей scRNA-seq клеток CSF доноров с MS и котроль; и по 5 повторностей scRNA-seq клеток PBMC доноров c MS и контроль
Препроцессинг каждой повторности:

Матрицы экспрессии фильтровались по минимальному колличеству клеток(3); минимальному и максимальному колличеству UMI (200, 2500); процент генов белков гемоглобина(10%); процент митхондриальных генов (10%)

Нормализация данных каждой повторности с помощью отрицательного биноминального распределения, реализованного во втроенном в seurat методе SCTransform

Затем повторности были интегрированны в один датасет методом IntegrateData. Получены 4 датасета: scf_ms, scf_cont, pbmc_ms, pbmc_cont




Далее была произведена кластеризация UMAP c первичным уменьшением размерости PCA
Результаты:

CSF_MS и CSF_cont:


![CSF_ms](https://github.com/Barabaika/Single_cell_csf/blob/main/csf_ms1.png)
![CSF_cont](https://github.com/Barabaika/Single_cell_csf/blob/main/csf_cont1.png)


PBMC_MS и PBMC_cont:


![PBMC_ms](https://github.com/Barabaika/Single_cell_csf/blob/main/pbmc_ms2.png)
![PBMC_cont](https://github.com/Barabaika/Single_cell_csf/blob/main/pbmc_cont3.png)

Найдены маркеры для каждого кластера 
