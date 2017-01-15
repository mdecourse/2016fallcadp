#coding: utf-8
# STL viewer 原始檔案來自
# University of Wuppertal - http://mbi-wiki.uni-wuppertal.de/wordpress/
# Modified by Uli Eggersmann
# Binary STL 資料讀取原始作者 Oliver Marks - http://www.linux.com
# 原始檔案僅讀取 Text STL 零件檔案
# 2011 Fall 由 KMOL 新增 Binary STL 零件檔案讀取

from visual import scene, color, materials, faces, points
import os, struct

#file ="ritzel.stl"

file ="binary.stl"

scene.width = 400
scene.height = 400
scene.background = color.white # black
# 視窗標題取自 cvisual.pyd, 不可使用中文
scene.title = "STLViewer in VPython"
 
print ("利用滑鼠右鍵旋轉")
print ("滑鼠左右鍵同時按下後移動, 可以縮放畫面")
 
# Read STL file, only use vertex-line with xyz coordinates
list = []

#load stl file detects if the file is a text file or binary file
def load_stl(filename):
    #read start of file to determine if its a binay stl file or a ascii stl file
    fp=open(filename,'rb')
    header=fp.read(80)
    filetype=header[0:5]
    # 這裡必須要能夠分辨二位元字串與文字字串
    #print (type(filetype))
    #print (filetype)
    fp.close()
    
    # for Python 3
    if filetype==b'solid':
    # for Python 2
    #if filetype=='solid':
        print ("讀取文字檔案格式:"+str(filename))
        load_text_stl(filename)
    else:
        print ("讀取二位元檔案格式:"+str(filename,))
        load_binary_stl(filename)

#load binary stl file check wikipedia for the binary layout of the file
#we use the struct library to read in and convert binary data into a format we can use
def load_binary_stl(filename):
    '''
    二位元 STL 檔案格式如下:
    檔案標頭共有 80 個字元(bytes), 內容通常省略, 但是內容不可使用 solid, 以免與文字檔案 STL 混淆
    UINT8[80] – Header
    UINT32 – Number of triangles (I:佔 4 bytes 的 unsigned integer)

    foreach triangle
    REAL32[3] – Normal vector (f:每一座標分量為一佔 4 bytes 的 float, 共佔 12 bytes)
    REAL32[3] – Vertex 1
    REAL32[3] – Vertex 2
    REAL32[3] – Vertex 3
    UINT16 – Attribute byte count (H:兩個 bytes 的 unsigned short, 表示 attribute byte count)
    end

    '''
    global list
    
    fp=open(filename,'rb')
    header=fp.read(80)

    triangle_number = struct.unpack('I',fp.read(4))[0]
    count=0
    while True:
        try:
            p=fp.read(12)
            if len(p)==12:
                n=[struct.unpack('f',p[0:4])[0],struct.unpack('f',p[4:8])[0],struct.unpack('f',p[8:12])[0]]
              
            p=fp.read(12)
            if len(p)==12:
                p1=[struct.unpack('f',p[0:4])[0],struct.unpack('f',p[4:8])[0],struct.unpack('f',p[8:12])[0]]
                list.append(p1)
            p=fp.read(12)
            if len(p)==12:
                p2=[struct.unpack('f',p[0:4])[0],struct.unpack('f',p[4:8])[0],struct.unpack('f',p[8:12])[0]]
                list.append(p2)
            p=fp.read(12)
            if len(p)==12:
                p3=[struct.unpack('f',p[0:4])[0],struct.unpack('f',p[4:8])[0],struct.unpack('f',p[8:12])[0]]
                list.append(p3)
            # 使用 count 來計算三角形平面個數
            # triangle_number 為 STL 檔案中的三角形個數
            count += 1
            # 在前面 12*4 個 bytes 的 normal 與三個點資料後, 為
            # 一個 2 bytes 長的 unsigned short, 其值為零, 為 attribute
            fp.read(2)
            # 讀完所有三角平面後, 即跳出 while
            if count > triangle_number:
                break
        except EOFError:
            break
    fp.close()


def load_text_stl(filename):
    global list
    for dataline in open(filename,"r").readlines():
        if not dataline.strip(): # skip blank lines
            continue
        field = dataline.split() # split with no argument makes the right place!
        if field[0] == "vertex":
            list.append([float(x) for x in field[1:4]])
            #print (list)
            #break
            #for x in field[1:4]:
                #print(x)

load_stl(os.path.abspath('')+'/'+file)


# Graphics
model = faces(pos=list, color=(0.8,0.8,0.8),
              material=materials.plastic) # creates triangles
# 請注意, 這裡並沒有使用 STL 檔案中的平面 normal, 而是利用 VPython make_normals() 產生
model.make_normals() # creates plane normals
model.smooth(0.93) # smooths the edges
# = AllepunkteSTL points (pos = list, size = 3, color = Color.Black) # generates points 