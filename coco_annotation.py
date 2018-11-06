import json
from collections import defaultdict
import os
from os import getcwd

sets = ['train2017','train2014','val2017','val2014']

for imagesss in sets:
    wd = getcwd()
    name_box_id = defaultdict(list)
    id_name = dict()
    f = open(
        'COCO/annotations/instances_%s.json'%(imagesss)  ,
        encoding='utf-8')
    data = json.load(f)
    annotations = data['annotations']
    for ant in annotations:
        id = ant['image_id']
        name = '%s/COCO/%s/%012d.jpg' %(wd,imagesss,id) 
        cat = ant['category_id']
        if cat >= 1 and cat <= 11:
            cat = cat - 1
        elif cat >= 13 and cat <= 25:
            cat = cat - 2
        elif cat >= 27 and cat <= 28:
            cat = cat - 3
        elif cat >= 31 and cat <= 44:
            cat = cat - 5
        elif cat >= 46 and cat <= 65:
            cat = cat - 6
        elif cat == 67:
            cat = cat - 7
        elif cat == 70:
            cat = cat - 9
        elif cat >= 72 and cat <= 82:
            cat = cat - 10
        elif cat >= 84 and cat <= 90:
            cat = cat - 11
        #测试是否有错误类型  
        if cat>79 or cat<0:
            print('wrong!\n') 
        #行人是第一类，cat=0,将其输出 
        if cat==0:
            name_box_id[name].append([ant['bbox'], cat])
    f = open('coco_%s.txt'%(imagesss), 'w')
    for key in name_box_id.keys():
        f.write(key)
        box_infos = name_box_id[key]
        for info in box_infos:
            x_min = int(info[0][0])
            y_min = int(info[0][1])
            x_max = x_min + int(info[0][2])
            y_max = y_min + int(info[0][3])
            box_info = " %d,%d,%d,%d,%d" % (
                x_min, y_min, x_max, y_max, int(info[1]))
            f.write(box_info)
        f.write('\n')
    f.close()
os.system("cat coco_train2017.txt coco_train2014.txt coco_val2017.txt coco_val2014.txt > coco_person.txt")
os.system("rm -rf coco_train2017.txt coco_train2014.txt coco_val2017.txt coco_val2014.txt ")

