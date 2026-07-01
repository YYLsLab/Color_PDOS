#!/bin/bash
######## 输入参数 ###############
atom_per_layer=41  # 每层多少原子i
num_layer=10 # 总共分多少层，i*num_layer 应该略小于总原子数
plot_direction="x"  #沿着哪个轴画图 x/y/z

# 去除screen_direct方向上坐标大于(x,y,z)/小于(-x,-y,-z)screen_coor的原子PDOS
screen_direct="z"  # 空值不执行筛选
screen_coor=0.5
################################
if [[ -n "$screen_direct" ]];  then  # 是否过滤
  screen_sign=${screen_direct:0:1} # 确定大于/小于
  # 确认screen列
  if [[ $screen_direct =~ "x" ]]; then
    screen_col=2
  elif [[ $screen_direct =~ "y" ]]; then
    screen_col=3
  elif [[ $screen_direct =~ "z" ]]; then
    screen_col=4
  fi
fi
natom=$( sed -n '1p' atom.config |  awk  '{print $1}' )  # 从atom.config文件中读取原子总数
atom_line=$( echo " 6 + $natom " | bc ) # last line of atom
lastline=$( echo "$natom  + 6 + 1 " | bc ) # last +1 line of atom
sed -i "$lastline,\$d" atom.config  # 保留原子坐标信息，删除受力等其他信息
spin=$( grep -i " spin " REPORT |awk '{print $NF}' ) # 从REPORT文件中判定自旋开启情况
echo  "spin = $spin"
rm color_DOS_for_origin.dat
cat -n atom.config >atom.config_withorder   #添加顺序
# 确认排序列
if [[ $plot_direction == "x" ]]; then
  sort_wiorder=3
elif [[ $plot_direction == "y" ]]; then
  sort_wiorder=4
elif [[ $plot_direction == "z" ]]; then
  sort_wiorder=5
fi
sort_woorder=$( echo " $sort_wiorder - 1 " | bc )
vim -c ' silent   :7,$!sort -n -k '"$sort_wiorder"  -c ':wq' atom.config_withorder  #按照指定排序-带顺序
awk '{print $1}' atom.config_withorder>order   #输出排序后的顺序
cut -f2- atom.config_withorder >atom.config    #移除顺序


if [[ $spin == "2" ]]; then
  for j in spinup  spindown; do
    for i in `seq -w 1 1 $num_layer`; do  #总共分为多少层，i*layer应该略小于总原子数
      mkdir PDOS-$i
      partial_and_coor.x   $i  $atom_per_layer $plot_direction  #生成每一层的atom.config和平均坐标
      mv atom.config_1 atom.config_mixed
      
      paste order atom.config_mixed >atom.config_mixed_wioder #合并顺序原子
      vim -c ' silent   :7,$!sort -n -k 1' -c ':wq' atom.config_mixed_wioder #顺序还原
      cut -f2- atom.config_mixed_wioder >atom.config #删除顺序标号
      
      # 执行筛选 =================================
      if [[ -n "$screen_direct" ]];  then
        for k in `seq 7 $atom_line` ;do
          coor_val=$( sed -n "${k}p" atom.config | awk -v c="$screen_col" '{print $c}')
          if [[ $screen_sign == "-" ]];then
            if [[ $(echo "$coor_val < $screen_coor" | bc) -eq 1 ]]; then   sed -i "${k}s/1$/0/" atom.config ; fi
          else
            if [[ $(echo "$coor_val > $screen_coor" | bc) -eq 1 ]]; then   sed -i "${k}s/1$/0/" atom.config ; fi
          fi
        done
      fi
      # ==============================================================
      plot_DOS_interp.x  > /dev/null  # pwmat的DOS脚本，生成DOS.totalspin/spinup/spindown等文件\
      cp DOS.$j DOS.totalspin
      cp atom.config DOS.* PDOS-$i
      
      vim -c ' silent   :7,$!sort -n -k '"$sort_woorder"  -c ':wq' atom.config  #不带顺序，已排序的原子-准备下次循环
      
      merge_DOS.x  $i   # 读取DOS.totalspin，生成merge_DOS.dat
    done
    merge_all_origin.x    # 合并merge_DOS.dat数据，生成merge_all_origin.dat数据
    mv merge_all_origin.dat merge_all_origin_${j}_${plot_direction}${screen_direct}.dat
  done

else
  for i in `seq -w 1 1 $num_layer`; do  #总共分为多少层，i*layer应该略小于总原子数
    mkdir PDOS-$i
    partial_and_coor.x   $i  $atom_per_layer $plot_direction    #生成每一层的atom.config和平均坐标
    mv atom.config_1 atom.config_mixed
    
    paste order atom.config_mixed >atom.config_mixed_wioder #合并顺序原子
    vim -c ' silent   :7,$!sort -n -k 1' -c ':wq' atom.config_mixed_wioder #顺序还原
    cut -f2- atom.config_mixed_wioder >atom.config #删除顺序标号
    # 执行筛选 =================================
    if [[ -n "$screen_direct" ]];  then
      for k in `seq 7 $atom_line` ;do
        coor_val=$( sed -n "${k}p" atom.config | awk -v c="$screen_col" '{print $c}')
        if [[ $screen_sign == "-" ]];then
          if [[ $(echo "$coor_val < $screen_coor" | bc) -eq 1 ]]; then   sed -i "${k}s/1$/0/" atom.config ; fi
        else
          if [[ $(echo "$coor_val > $screen_coor" | bc) -eq 1 ]]; then   sed -i "${k}s/1$/0/" atom.config ; fi
        fi
      done
    fi
    # ==============================================================
    plot_DOS_interp.x > /dev/null  # pwmat的DOS脚本，生成DOS.totalspin/spinup/spindown等文件\
    cp atom.config DOS.* PDOS-$i
    
    vim -c ' silent   :7,$!sort -n -k '"$sort_woorder"  -c ':wq' atom.config  #不带顺序，已排序的原子-准备下次循环
    
    merge_DOS.x  $i   # 读取DOS.totalspin，生成merge_DOS.dat
  done
  merge_all_origin.x    # 读取merge_DOS.dat，生成merge_all_origin.dat数据
  mv merge_all_origin.dat merge_all_origin_${plot_direction}${screen_direct}.dat
fi

paste order atom.config > atom.config_withorder #合并顺序-原子
vim -c ' silent      :7,$!sort -n -k 1' -c ':wq' atom.config_withorder #顺序还原
cut -f2- atom.config_withorder >atom.config #删除顺序标号
