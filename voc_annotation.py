import xml.etree.ElementTree as ET
import os
from os import getcwd

sets=[('2012', 'train'), ('2012', 'val'),('2007', 'train'), ('2007', 'val'), ('2007', 'test')]
classes = ["person"]

wd = getcwd()

for year, image_set in sets:
    image_ids = open('VOCdevkit/VOC%s/ImageSets/Main/%s.txt'%(year, image_set)).read().strip().split()
    list_file = open('%s_%s.txt'%(year, image_set), 'w')
    for image_id in image_ids: 
        in_file = open('VOCdevkit/VOC%s/Annotations/%s.xml'%(year, image_id))
        tree=ET.parse(in_file)
        root = tree.getroot()

        is_exit = 0
        for obj in root.iter('object'):
            difficult = obj.find('difficult').text
            cls = obj.find('name').text
            if cls in classes : 
                list_file.write('%s/VOCdevkit/VOC%s/JPEGImages/%s.jpg'%(wd, year, image_id))    
                is_exit = 1            
                break 

        if(is_exit): 
            for obj in root.iter('object'):
                difficult = obj.find('difficult').text
                cls = obj.find('name').text
                if cls not in classes or int(difficult)==1: 
                    continue
                cls_id = classes.index(cls)
                xmlbox = obj.find('bndbox')
                b = (int(xmlbox.find('xmin').text), int(xmlbox.find('ymin').text), int(xmlbox.find('xmax').text), int(xmlbox.find('ymax').text))
                list_file.write(" " + ",".join([str(a) for a in b]) + ',' + str(cls_id))
            list_file.write('\n')

list_file.close()


os.system("cat 2007_train.txt 2007_val.txt 2007_test.txt 2012_train.txt 2012_val.txt > voc_person.txt")
os.system("rm -rf 2007_train.txt 2007_val.txt 2007_test.txt 2012_train.txt 2012_val.txt")

