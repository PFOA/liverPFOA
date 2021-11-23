export R_LIBS="/szrmyy/wangjgLab/scRNA/chenjh/platform/software/lib/R/library:/szrmyy/wangjgLab/scRNA/chenjh/platform/backup/R/R-4.1.0/library:$R_LIBS" 
Rscript   extract_epi_seurat.R  
  /szrmyy/wangjgLab/scRNA/chenjh/platform/software/project/cellrank/bin/python  make_anndata.py  
 /szrmyy/wangjgLab/scRNA/chenjh/platform/software/project/cellrank/bin/python    create_caches.py 
