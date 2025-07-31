python ../design/barcode_design.py \
    --barcodeLen 20 \
    --polymerLen 3 \
    --minGCRate 0.4 \
    --maxGCRate 0.6 \
    --dictRepeats "GC:3,CG:3,AT:3,TA:3" \
    --innerPairLen 4 \
    --outerPairLen 8 \
    --minEditDistance 9 \
    --numCandidateBarcode 500 \
    --numGeneratedBarcode 96 \
    --numProcess 16 \
    --outputSeqFile barcodes.fa \
    --outputDistanceFile distances.txt \
    --logFile run.log


python ../design/barcode_library.py \
    barcodes.fa library.txt 50
