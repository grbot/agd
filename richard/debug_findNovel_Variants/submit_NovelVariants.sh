
#!/bin/bash
#PBS -N Novel-Variants-Calc

python /spaces/gapw/diversity/richard/findNovelVariants_Benin.py -p /spaces/gapw/diversity/mamana/VCF_POP/ -o /spaces/gapw/diversity/richard/Novel_variants/Benin_Novel_variants.txt

python /spaces/gapw/diversity/richard/findNovelVariants_Botswana.py -p /spaces/gapw/diversity/mamana/VCF_POP/ -o /spaces/gapw/diversity/richard/Novel_variants/ Botswana_Novel_variants.txt

python /spaces/gapw/diversity/richard/findNovelVariants_Burkina.py -p /spaces/gapw/diversity/mamana/VCF_POP/ -o /spaces/gapw/diversity/richard/Novel_variants/ Burkina_Novel_variants.txt

python /spaces/gapw/diversity/richard/findNovelVariants_Cameroon.py -p /spaces/gapw/diversity/mamana/VCF_POP/ -o /spaces/gapw/diversity/richard/Novel_variants/Cameroon_Novel_variants.txt

python /spaces/gapw/diversity/richard/findNovelVariants_Ghana.py -p /spaces/gapw/diversity/mamana/VCF_POP/ -o /spaces/gapw/diversity/richard/Novel_variants/ Ghana_Novel_variants.txt

python /spaces/gapw/diversity/richard/findNovelVariants_Mali.py -p /spaces/gapw/diversity/mamana/VCF_POP/ -o /spaces/gapw/diversity/richard/Novel_variants/Mali_Novel_variants.txt

python /spaces/gapw/diversity/richard/findNovelVariants_Nigeria.py -p /spaces/gapw/diversity/mamana/VCF_POP/ -o /spaces/gapw/diversity/richard/Novel_variants/Nigeria_Novel_variants.txt

python /spaces/gapw/diversity/richard/findNovelVariants_Zambia.py -p /spaces/gapw/diversity/mamana/VCF_POP/ -o /spaces/gapw/diversity/richard/Novel_variants/Zambia_Novel_variants.txt
