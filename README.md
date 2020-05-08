# 示例说明

该模型代码主要来自于GBEHM模型以及相关的陆面过程模型CoLM，并耦合了风吹雪模型PBSM。在WSP文件中给出一个模拟示例，区域为新疆玛纳斯河流域。
以下为相关说明事项。文字来自于在兰州寒旱所举行的模型培训班笔记，主要撰写者为吴雪娇，并由马佳培、雷华锦整理。

### 环境安装与配置

a. ArcGIS 10.5 安装(流域提取需使用 python27 的 f2py) 注意事项： 1）保证电脑是 window10 系统，且没有系统文件被损坏。 2）如果原电脑有 10.5 以下版本需卸载，则一定要完全卸载，才可安装 10.5. 3）在 python27 目录下用 pip 更新 numpy,并安装 tqdm。 4）环境变量配置。 环境变量：C:\Python27\ArcGIS10.5 系统变量：C:\Python27\ArcGIS10.5\python.exe;C:\Python27

b. MinGW 安装(为了让 python 调用 fortran) 注意事项：

1) 安装时可能需要VPN，才能安装成功。

2) 环境变量：C:\MinGW\bin

c. 在 windows 安装 Linux 环境 注意事项：

1) 设置->更新和安全->开发者选项，选择开发模式。

2) 控制面板->程序->打开关闭 windows 窗口，选择 windows 子系统 Linux 选项。

3) Ubuntu，选择并安装。设置用户名和密码。

4) 在 linux 系统下更新此系统。注意，如果更新不成功可能是源码库调用 不成功。就要换源码库。 5) 安装 cdo, nco, dos2unix, gfortran, python3-dev 备注：另建议安装软件有， Notepad++, PyCharm professional, PyNcView, Panoply, QGIS.

### 数据准备

a. 流域数据提取(Windows 环境下运行) 1）先编译 .\compile_Fhorton.bat 2）后运行 python SubBaisnExtract_CMD.py 注意事项： 1）流域提取需进行两次，第一次运行使用必选关键字，-d, 100 米分辨率的 DEM；-p, 流域出山口的控制点 shp 文件位置，-w 输出文件路径。根据结 果制作 subcatchment.dat 文件。第二次运行，除上述关键字外，还需添加 s, subchatchment.dat 文件路径及 –l, landuse.tif 文件路径。注意：土地利用 数据需要 USGS 分类方法的数据。

b. 气象数据准备（Linux 环境下运行） 气象数据需要两类，一类是分析数据（空气温度，水平和垂直风速，相对 湿度，气压，比湿）；另一类是预报数据（降水，长短辐射）。1h 的时间分 辨率和 1km 的空间分辨率。每月一个文体。需要注意的是要严格按照模型 要求的单位来准备所有数据。数据是 nc 格式的，建议用 cdo 来处理。 ERA_Interim 数据可用./ERA_convert.sh 来处理，进行单位转换，会产生 ERA_forecat.nc 和 ERA_reanalysis.nc 两个文件；再用。./resample_metro.sh 时 空插值，这里注意，需要流域提取出来的 coord.nc 这个文件，并用 cdo griddes 转换成 coord.txt。需要注意的是，ERA 资料由于启始时间问题，如处理 2005 年 1-12 月的数据，需准备 2004.12-2006.1 这 14 个月的文件。

c. 土壤数据提取 示例是南京土壤所制作的全国土壤分布图。./resample_soil.sh 这里也需要 coord.txt 文件。

### 参数率定 
有 6 个参数需要率定，根据不同流域的实际情况。包括：表层土壤深度；承压 水层深度（潜水层+含水层）；土壤孔隙度；含水层储水量；蒸散调节系数； 地下水传导率。这些都在主程序中。

### 模型运行 
1）在 code 中，编译./linux_f2py_bat, 来生产 pygbhm.so 和 pyice.so 两个文件， 如果是 python3 则生成的是.so 和.so 两个文件。要把这两个文件复制到主程序 main.py 所在的文件夹。 2）主程序运行 python3 main.py 这里可能需要装 python3-numpy; python3-pandas; python3-netcdf; gdal-bin.

在主程序中修改气象输入和模拟结果路径；在 AuxTools.py 中修改气象数据名 称。 GDAL 是个非常强大的工具，ArcGIS 的底层也是调用此工具。 注意：模型开始运行的时候，第一年需要运行 4 次预热。第五次才开始得到 第一年的结果。
