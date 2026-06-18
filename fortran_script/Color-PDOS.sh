#!/bin/bash
rm color_DOS_for_origin.dat
cat -n atom.config >atom.config_withorder
vim -c ':7,$!sort -n -k 5' -c ':wq' atom.config_withorder
awk '{print $1}' atom.config_withorder>order
cut -f2- atom.config_withorder >atom.config
atom_per_layer=25  #每层多少原子
for i in `seq -w 1 1 10`  #总共分为多少层，i*layer应该略小于总原子数
do
  mkdir PDOS-$i
  ./partial_and_coor.x   $i  $atom_per_layer   #生成每一层的atom.config和平均坐标zsffdd 
  mv atom.config_1 atom.config_mixed
  
  paste order atom.config_mixed >atom.config_mixed_wioder #合并顺序原子
  vim -c ':7,$!sort -n -k 1' -c ':wq' atom.config_mixed_wioder #顺序还原
  cut -f2- atom.config_mixed_wioder >atom.config #删除顺序标号
  
  plot_DOS_interp.x
  cp atom.config DOS.* PDOS-$i
  
  vim -c ':7,$!sort -n -k 4' -c ':wq' atom.config #再次打乱 准备下次循环
  ./merge_DOS.x  $i  
done
./merge_all_origin.x  
paste order atom.config > atom.config_withorder #合并顺序-原子
vim -c ':7,$!sort -n -k 1' -c ':wq' atom.config_withorder #顺序还原
mv atom.config_withorder atom.config
