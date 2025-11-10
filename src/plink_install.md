# installing PLINK
I am curious about using PLINK to filter for long-range LD. Since it was not present in `/projects/gatins/programs_explorer/`, I installed it.

new folder

`mkdir -p plink
cd plink`

then install latest version
`wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip`

unzip
`unzip plink_linux_x86_64_20231211.zip`

make executeable 
`chmod +x plink`

add to session path if you'd like!
`export PATH=$PATH:/projects/gatins/programs_explorer/plink`
