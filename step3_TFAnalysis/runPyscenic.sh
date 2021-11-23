#/usr/bin/Rscript    getsampleCSV.R
export PYTHONPATH="/home/chenjunhui/anaconda3/envs/pyscenic/lib/python3.6/site-packages:/home/chenjunhui/miniconda3/lib/python3.6/site-packages"
/home/chenjunhui/anaconda3/envs/pyscenic/bin/python3  produceLoom.py
/home/chenjunhui/anaconda3/envs/pyscenic/bin/pyscenic     grn    sample.loom   /szrmyy/wangjgLab/scRNA/chenjh/P0003_ratLung/08.TF/SCENICprotocol/example/allTFs_mm.txt   -o adj.csv --num_workers 20
/home/chenjunhui/anaconda3/envs/pyscenic/bin/pyscenic   ctx  adj.csv   /public/group/chenjunhui/scRNA/TF/pyscienicDatabase/humanFeature/cisTarget/mm9-500bp-upstream-7species.mc9nr.feather   /public/group/chenjunhui/scRNA/TF/pyscienicDatabase/humanFeature/cisTarget/mm9-tss-centered-10kb-7species.mc9nr.feather  --annotations_fname   /public/group/chenjunhui/scRNA/TF/pyscienicDatabase/motifs-v9-nr.mgi-m0.001-o0.0.tbl      --expression_mtx_fname  sample.loom  --output reg.csv    --mask_dropouts   --num_workers 20

/home/chenjunhui/anaconda3/envs/pyscenic/bin/pyscenic  aucell \
   sample.loom   \
    reg.csv \
    --output  Samplef_loom_path_scenic_output.loom \
    --num_workers 20