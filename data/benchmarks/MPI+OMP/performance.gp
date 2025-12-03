reset
set terminal pngcairo size 1200,600 enhanced font 'Arial,12'
set output 'performance.png'

# 设置柱状图样式
set style data histograms
set style histogram rowstacked
#set style histogram clustered gap 1
set style fill solid border -1
set boxwidth 0.6

# 设置坐标轴
set title "Pure MPI(size 2000 * 1000)"
set xlabel "Number of Nodes"
set ylabel "Time (s)"
set grid ytics

# 使用不同颜色
set key outside

# 绘图
plot 'time.txt' using 3:xtic(1) title "T_init", \
     '' using 4 title "T_compute", \
     '' using 5 title "T_communicate", \
     '' using 6 title "T_output"

# 第二张图：Efficiency & Speedup
set output 'efficiency_speedup.png'
set terminal pngcairo size 1200,600 enhanced font 'Arial,12'

# 计算 speedup 和 efficiency
# speedup = T_total(1 node) / T_total(N)
# efficiency = speedup / N
T1 = 3086.582649028  # 2节点的T_total，如果想以1节点为基准需调整数据

set title "Efficiency and speedup(size 2000 * 1000)"
set xlabel "Number of Nodes"
set ylabel "Speedup / Efficiency"
set grid ytics
set ytics nomirror
set y2label "Efficiency"
set y2tics

set key outside

plot 'time.txt' using ( $1 ):(T1/$2) with linespoints lt 1 lw 2 pt 7 title "Speedup", \
     '' using ( $1 ):(T1/$2/$1) axes x1y2 with linespoints lt 2 lw 2 pt 5 title "Efficiency"
