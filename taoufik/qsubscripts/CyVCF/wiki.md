# Wiki
Script should be executed under the cluster.

  - CyVCF need to be installed
  - Run script with full paths for both VCF and result folders 

### CyVCF installation

```sh
$ sudo pip install cyvcf2
$ sudo easy_install cyvcf2
```
### Running the script
```sh
$ ./NovelVariantsAnnotation.sh /Path/To/VCF/Folder/ /Path/To/Results/Folder
$ sh NovelVariantsAnnotation.sh /Path/To/VCF/Folder/ /Path/To/Results/Folder
```
### Jobs infromations
For a file  VCF-file.vcf.gz :

Job name :  
```sh 
Job-Novel-Variants-in-VCF-file.vcf.gz
```
Job error output :  
```sh 
VCF-file.vcf-error.log
```
Job stdout output :  
```sh 
VCF-file.vcf.gz-out.log
```

