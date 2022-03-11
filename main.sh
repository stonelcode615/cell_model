Ion_index=0           #   0-Na; 1-K; 2-Rb; 3-Cs
opt=1                 # 0 compute; 1 analysis

home_folder=`pwd`

Ion_names=(na k rb cs)
Ion_sizes=(1.33 1.77 1.92 2.14)
Ion_gauss=(0.67 0.96 1.12 1.47)
Ion_Bis=(0.398 -0.340 0.000 -0.790)

Ion_g1s=(47.59 10.44 0.00 13.53)

paths=(0-Na 1-K 2-Rb 3-Cs)

na_nf=(12 14 16 18)
k_nf=(12 14 16 18)
rb_nf=(14 22)
cs_nf=(22 )

na_sv=(1.002 1.022 1.041 1.055)
k_sv=(0.973 0.996 1.015 1.031)
rb_sv=(0.866 0.952)
cs_sv=(0.873 )

na_H0=(0.0383 0.0356 0.0335 0.0320)
k_H0=(0.0282 0.0257 0.0229 0.0206)

name=${Ion_names[$Ion_index]}
size=${Ion_sizes[$Ion_index]}
gauss=${Ion_gauss[$Ion_index]}
Bi=${Ion_Bis[$Ion_index]}
g1=${Ion_g1s[$Ion_index]}

echo $name $size $gauss $Bi $g1

NUM=$(eval echo "\${#${name}_nf[@]}")
echo $NUM

for n in `seq 0 1 $[NUM-1]`
do
    sv=$(eval echo \${${name}_sv[$n]})
    nf=$(eval echo \${${name}_nf[$n]})
    nc=`echo "$nf-1" | bc`

    mu10=0.000

    mu20=0.000
    g0=0.000
    H0=$(eval echo \${${name}_H0[$n]})

    mu30=0.000

    work_folder="./$nf"
    if [ $opt == 0 ]; then
        if [ ! -d $folder_path ]; then
            mkdir $folder_path
        fi
        echo $n $nf $sv
        cp -r ./package $work_folder
        cd $work_folder
        sed -i -e "s/@name/$name/g" -e "s/@nf/$nf/g"                                        cell.job
        sed -i -e "s/@name/$name/g" -e "s/@nc/$nc/g" -e "s/@sv/$sv/g"                       setup.py
        sed -i -e "s/@size/$size/g" -e "s/@gauss/$gauss/"                                   setup.py
        sed -i -e "s/@mu10/$mu10/g" -e "s/@g1/$g1/g" -e "s/@Bi/$Bi/g"                       compute.py
        sed -i -e "s/@mu20/$mu20/g" -e "s/@g0/$g0/g" -e "s/@H0/$H0/g"                       compute.py
        sed -i -e "s/@mu30/$mu30/g"                                                         compute.py
        qsub cell.job
        cd $home_folder
    fi
    if [ $opt == 1 ]; then
        cd $work_folder
        python3 0log-generate.py
    fi
    cd $home_folder
done








