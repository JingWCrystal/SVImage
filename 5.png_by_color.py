#coding:utf-8 
from PIL import Image
import os
def cluster_dic(img_path,image,res_path,res_file):
	# 把带有breakpoints的小图片拿出来分析，每行，有颜色的像素标为1，没颜色的标为0
	img=Image.open(img_path+image)
	img_width=img.size[0]
	img_height=img.size[1]
	newIm = Image.new ("RGB", (img_width,img_height),(255,255,255))
	# 把图片的每行变成0,1,2,3
	img_arr=[]
	for line in range(img_height):
		tag=0 # 优先级： clip 3 > discorcondant 2 > corcondant 1 > no base 0
		for wid in range(img_width):
			pixel=img.getpixel((wid,line))
			if pixel[0]==255 and not pixel[1]==255 and not pixel[2]==255: #含有clipped 的优先，
				tag=3
			elif not pixel[0]==255 and not pixel[1]==255 and pixel[2]==255: #discorcondant 的
				if tag <3 :
					tag=2
			elif pixel[0]==255 and not pixel[1]==255 and not pixel[2]==255: #corcondant 的
				if tag <2 :
					tag=1
		img_arr.append(tag)
	list_rearrange=[]
	list_clip=[]
	list_discorcondant=[]
	list_corcondant=[]
	list_AllSpace=[]
	for i in range(len(img_arr)):
		if img_arr[i]==0:
			list_AllSpace.append(i) # 行号
		elif img_arr[i]==1:
			list_corcondant.append(i) # 行号
		elif img_arr[i]==2:
			list_discorcondant.append(i) # 行号
		elif img_arr[i]==3:
			list_clip.append(i) # 行号
	list_rearrange.extend(list_clip)
	list_rearrange.extend(list_discorcondant)
	list_rearrange.extend(list_corcondant)
	list_rearrange.extend(list_AllSpace)
	rm_blank=False
	for row in range(img_height):
		for col in range(img_width):
			red,green,blue=img.getpixel((col,list_rearrange[row]))[0:3]
			if not (red==255 and green==255 and blue==255):
				rm_blank=True
			newIm.putpixel((col,row),(red,green,blue))
	if rm_blank:
		newIm.save(res_path+image,"PNG")

	imf=open(res_file,'a')
	imf.write(image+"\t"+'0'+"\n") #[3] 不区分纯合杂合，[2]是区分
	imf.close()


# waiting---------------------------------
def main():
	img_path='/mnt/hde/gao/wj/keras/left_0/'
	res_path='/mnt/hde/gao/wj/keras/left_bycolor/'
	res_file='/mnt/hde/gao/wj/keras/colored_left.txt'
	img_width=200
	img_height=100
	for file_tmp in os.listdir(img_path):
		if 'left' in str(file_tmp) and 'png' in str(file_tmp):
			# print file_tmp
			cluster_dic(img_path,file_tmp,res_path,res_file)

	img_path='/mnt/hde/gao/wj/keras/right_0/'
	res_path='/mnt/hde/gao/wj/keras/right_bycolor/'
	res_file='/mnt/hde/gao/wj/keras/colored_right.txt'
	img_width=200
	img_height=100
	for file_tmp in os.listdir(img_path):
		if 'right' in str(file_tmp) and 'png' in str(file_tmp):
			# print file_tmp
			cluster_dic(img_path,file_tmp,res_path,res_file)


if __name__ == '__main__':
    main()