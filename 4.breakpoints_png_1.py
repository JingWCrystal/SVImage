#coding:utf-8 

import pysam
from PIL import Image

def rearrange_string(read):
	bases=read.query_sequence   #原read序列的碱基们
	# print read.cigarstring
	# print bases
	new_bases=''               #重排后的序列
	new_base_quality=[]
	
	#1. M部分的照搬AGCT不变,颜色4种；new_bases插入原始agct,移动bases_checked
	#2. I部分的base忽视掉，不向new_bases中填充，移动bases_checked
	#3. D部分的统一成"d",空白，向new_bases插入‘d’,bases_checked不移动
	#4. 将S部分的碱基们统一成"s",空白，new_bases插入‘s’,移动bases_checked
	#   bases_checked 只在base上作为指针移动，所以只对S,I,M的碱基移动，D区域就不变 
	#   base_quality  只在有base的地方存在，移动指针同bases_checked。               
	bases_checked=0   
	read_len_count=0
	for cigar_portion in read.cigartuples:# cigar中只有0,1,2,4  M.I.D.S 
		if cigar_portion[0]==0:#M
			cigar_base = bases[bases_checked:(cigar_portion[1]+bases_checked)]
			new_bases=new_bases+cigar_base
			for M_num in range(cigar_portion[1]):
				#print read_len_count
				new_base_quality.append( min(read.query_alignment_qualities[read_len_count],read.query_qualities[read_len_count]) )
				read_len_count=read_len_count+1
			bases_checked=bases_checked+cigar_portion[1] 

		elif cigar_portion[0]==1:#I
			bases_checked=bases_checked+cigar_portion[1] 
			#print read_len_count
			for I_num in range(cigar_portion[1]):
				read_len_count=read_len_count+1
				#print read_len_count
		elif cigar_portion[0]==2:#D  注意，deletion往原有序列里加了空白base，这个时候有可能使read超出右边界
			cigar_base=''
			for  i in range(cigar_portion[1]):
				cigar_base=cigar_base+'d'	
				new_base_quality.append(0)
			new_bases=new_bases+cigar_base
		elif cigar_portion[0]==4 : 
		#S  注意，往原有序列里加了空白base，这个时候有可能使read超出右边界
			cigar_base=''
			for  i in range(cigar_portion[1]):
				cigar_base=cigar_base+'s'
				new_base_quality.append(-1)	
				#print read_len_count
				# read_len_count=read_len_count+1
			new_bases=new_bases+cigar_base
			bases_checked=bases_checked+cigar_portion[1]
		elif cigar_portion[0]==5 : # hard clip 是不出现在read base当中的
			cigar_base=''
			for  i in range(cigar_portion[1]):
				cigar_base=cigar_base+'s'
				new_base_quality.append(-1)	
				#print read_len_count
				# read_len_count=read_len_count+1
			new_bases=new_bases+cigar_base
			# bases_checked=bases_checked+cigar_portion[1]
	# print read.cigarstring
	# print bases
	# print new_bases
	# print 'len(new_bases) = ',
	# print len(new_bases)
	# print 'len(new_base_quality) = ',
	# print len(new_base_quality)
	return new_bases,new_base_quality

def read_can_shown(read,scan_l_pos,scan_r_pos):
	read_pos1=read.reference_start
	read_pos2=read.reference_start+read_infered_len(read)
	# print ('if (read_pos2 > scan_l_pos) and (read_pos1 > scan_r_pos)')
	# print read_pos2
	# print scan_l_pos
	# print read_pos1
	# print read_infered_len(read)
	if (read_pos2 > scan_l_pos) and (read_pos1 < scan_r_pos):
		res=True
		for cigar_portion in read.cigartuples:
			if not ((cigar_portion[0]==0) or (cigar_portion[0]==1) or (cigar_portion[0]==2) or (cigar_portion[0]==4) or (cigar_portion[0]==5)):
				res = False 
		return res
	else:
		return False

def read_corner_shown(read,scan_l_pos,scan_r_pos,new_bases):#边界的read可以显示
	read_pos1=read.reference_start
	read_pos2=read.reference_start+read_infered_len(read)
	if (read_pos1 < scan_l_pos) and (read_pos2 > scan_l_pos):
		if ('A' or 'G' or 'C' or 'T' or 'a' or 'g' or 'c' or 't') in new_bases[(scan_l_pos-read_pos1):len(new_bases)] :
			return True
		else:
			return False
	if (read_pos1 < scan_r_pos) and (read_pos2 > scan_r_pos):
		if ('A' or 'G' or 'C' or 'T' or 'a' or 'g' or 'c' or 't') in new_bases[0:(scan_r_pos-read_pos1)] :
			return True
		else:
			return False
	if (read_pos1 >= scan_l_pos) and (read_pos2 <= scan_r_pos):
		return True


def read_infered_len(read): # 本来应是cigar0,1,4的相加，改成M,D,S的相加，0,2,4
	infer_len=0
	for cigar_portion in read.cigartuples:
		if ((cigar_portion[0]==0) or (cigar_portion[0]==2) or (cigar_portion[0]==4)):
			infer_len=infer_len+cigar_portion[1]
	return infer_len

def is_empty(read_list):
	tag=True
	for li in read_list:
		if li!=[]:
			tag=False
			break
	return tag

def get_shortest_tail_row(read_list,scan_r_pos):
	if is_empty(read_list):
		return 0
	else:
		tail=scan_r_pos
		short_row=0
		for i in range(len(read_list)):
			if read_list[i]!=[]:
				if tail > read_list[i][-1][1] :
					tail= read_list[i][-1][1]
					short_row=i
		return short_row

def find_next_empty_row(read_list):
	row=0
	for i in range(len(read_list)):
		if read_list[i]==[]:
			row=i
			break
	return row

def read_to_dictionary(read_package, scan_r_pos,height):
	#初始化read_list,有height个值，每个值是一个list,
	#i是行号，read_list[i]中存放该行的reads元组[(pos1,pos2,read,base,...),...]
	dictionary={}
	read_list=[0 for i in range(height)] 
	for i in range(height):
		read_list[i]=[]
	row_ptr=0 #dic_ptr 从0到height-1
	for base_and_read in read_package:
		if row_ptr<height:
			base=base_and_read[0]
			quality=base_and_read[1]
			read=base_and_read[2]
			
			#开头是S的，如果补全了S，start位点应该-len(s),因为此时的read.reference_start指的是从M开始的位点
			if read.cigartuples[0][0]==4 : 
				read_pos1=read.reference_start - read.cigartuples[0][1]
				read_pos2=read_pos1+read_infered_len(read)
			else:
				read_pos1=read.reference_start
				read_pos2=read_pos1+read_infered_len(read)
			if read.is_paired:  # read not paired，do not consider this kind of read
				is_concordant=read.is_proper_pair
				if 'S' in read.cigarstring:
					is_clipped=True
				else:
					is_clipped=False

				if is_empty(read_list):
					read_list[row_ptr].append( (read_pos1,read_pos2,base,quality,is_clipped,is_concordant) )
				else: 
					row_ptr = get_shortest_tail_row(read_list,scan_r_pos)
					if read_pos1 >= read_list[row_ptr][-1][1]: # 当前read的头在最短read tail之后，和这个放到同一行
						read_list[row_ptr].append( (read_pos1,read_pos2,base,quality,is_clipped,is_concordant) )
					else:                                           # 当前read的头在最短read tail之前，就放到最下一行
						row_ptr=find_next_empty_row(read_list)
						read_list[row_ptr].append( (read_pos1,read_pos2,base,quality,is_clipped,is_concordant))

	for i in range(height):
		dictionary[i]=read_list[i]
	return dictionary

def draw_pgn(which_bp,dic,width,height,scan_l_pos,scan_r_pos,img_name):
	newIm = Image.new ("RGB", (width,height),(255,255,255))
	for key in range(height):
		for read_tuple in dic[key]:
			read_pos1=read_tuple[0]
			read_pos2=read_tuple[1]
			base=read_tuple[2]
			quality=read_tuple[3]
			is_clipped=read_tuple[4]
			is_concordant=read_tuple[5]
			col=read_pos1-scan_l_pos
			index_in_read=0
			for i in range(len(base)):
				if col >=0  and col < width :
					row=key
					# print row,col
					red,green,blue=get_RGB(which_bp,base[index_in_read],quality[index_in_read],is_clipped,is_concordant)
					newIm.putpixel((col,row),(red,green,blue))
					index_in_read=index_in_read+1
					col=col+1
				elif col<0:
					index_in_read=index_in_read+1
					col=col+1
	newIm.save(img_name,"PNG")

def get_RGB(which_bp,base,quality,is_clipped,is_concordant):
	# print (which_bp)
	if is_clipped :#红色的是clip
		red=255
		green=0
		blue=0
		if quality>1:
			green=green+255-6*quality
			blue=blue+255-6*quality
			return red, green, blue
		elif quality== 0 or quality== -1 :
			return 255,255,255
		# elif quality== 0:  # del 的区域
		# 	return 255,255,255
		# elif quality==-1:   # soft clipped 和
		# 	if which_bp=='left':
		# 		return 255,255,255
		# 	elif which_bp=='right':
		# 		return 255, 255, 255
	elif is_concordant:#正常read pair是绿色
		red=0
		green=255
		blue=0
		if quality>1:
			red=red+255-6*quality
			blue=blue+255-6*quality
			return red, green, blue
		elif quality==0 or quality== -1:  # soft clipped 和 del 的区域
			return 255,255,255
	elif not is_concordant:#不正常read pair是蓝色
		red=0
		green=0
		blue=255
		if quality>1:
			red=red+255-6*quality
			green=green+255-6*quality
			return red, green, blue
		elif quality==0 or quality== -1:  # soft clipped 和 del 的区域
			return 255,255,255
	
def get_range(b1,b2):
	tmp=str(abs(b2-b1))
	high_bit=int(tmp[0])+1
	new_tmp=str(high_bit)
	for i in range(1,len(tmp)):
		new_tmp=new_tmp+'0'
	new_tmp=int(new_tmp)
	width= new_tmp+400
	height=150
	scan_l_pos= int(b1-(new_tmp-(b2-b1))/2-200)
	scan_r_pos=scan_l_pos+width

	# width= new_tmp+400
	# height=100
	# scan_l_pos= int(b1-(new_tmp-(b2-b1))/2-200)
	# scan_r_pos=scan_l_pos+width
	return width,height,scan_l_pos,scan_r_pos


# waiting------------------针对100bp以上的deletion的断点作图---------------
def main():
	# 纯合的
	vcf_path='/mnt/hde/gao/wj/simulate/testSimuResult1/bp.txt'
	# vcf_path='/mnt/hde/gao/wj/vcf/tools_made_vcf/NA12878_chr10.benchmark.vcf'
	bam_path='/mnt/hde/gao/wj/simulate/testSimuResult5/500_DEL_1000.bam'
	# 直接做 200*100的图，不截取了
	width=200
	height=100

	for line in open(vcf_path):
		print (line)
		line=line.strip('\n')
		tmp=line.split(' ')
		# chrom=tmp[0]      # chr1-chr7
		chrom='11:30000-36000000' 
		bk1=int(tmp[0])
		bk2=int(tmp[1])
		label='1'

		samfile = pysam.AlignmentFile(bam_path,"rb")
		#左断点~
		img_tmp='/mnt/hde/gao/wj/keras/left/'
		img_name=img_tmp+tmp[0]+'_'+tmp[1]+'.'+label+'.left.png'
		scan_l_pos=int(bk1-width/2)
		scan_r_pos=scan_l_pos+width
		# print (str(scan_l_pos))
		# print (str(scan_r_pos))
		read_package=[]
		c=0
#测试---------------------------------------------
		# img_tmp='/mnt/hde/gao/wj/SV_image/'
		# scan_l_pos=38816545
		# scan_r_pos=38818280
		# chrom='chr10'
		# label='1'
		# img_name='chr10.38816545-38816945.left.png'
		# width=1735
		# height=100
#测试---------------------------------------------		
		for read in samfile.fetch(chrom, scan_l_pos,scan_r_pos):
			if (read.cigarstring != None) and read_can_shown(read,scan_l_pos,scan_r_pos) and read.mapping_quality >10:
				new_bases,new_base_quality=rearrange_string(read)
				if read_corner_shown(read,scan_l_pos,scan_r_pos,new_bases):
					read_package.append((new_bases,new_base_quality,read))
					c=c+1
		dic=read_to_dictionary(read_package, scan_r_pos,height)
		print ('left: ok before drawn')
		if not c==0 : #图片上有read
			draw_pgn('left',dic,width,height,scan_l_pos,scan_r_pos,img_name)
			imf=open("/mnt/hde/gao/wj/keras/left/homo_left.txt",'a')
			imf.write(tmp[0]+'_'+tmp[1]+'.'+label+'.left.png'+"\t"+label+"\n")
			imf.close()

		# else : # 是一张空白图片，就不画这个图了！
		# 	imf=open("/mnt/hde/gao/wj/SV_image/l_blank_region.txt",'a')
		# 	imf.write(tmp[0]+"\t"+tmp[1]+"\t"+tmp[2]+"\t"+label+"\n")
		# 	imf.close()

		# # 右断点~
		img_tmp='/mnt/hde/gao/wj/keras/right/'
		img_name=img_tmp+tmp[0]+'_'+tmp[1]+'.'+label+'.right.png'
		scan_l_pos=int(bk2-width/2)
		scan_r_pos=scan_l_pos+width
# 测试----------------------------------------------
		# img_tmp='/mnt/hde/gao/wj/SV_image/'
		# scan_l_pos=38818080
		# scan_r_pos=38818480
		# chrom='chr10'
		# label='1'
		# img_name='chr10.38816545-38816945.right.png'
		# width=400
		# height=100
# 测试----------------------------------------------
		# print (str(scan_l_pos))
		# print (str(scan_r_pos))
		read_package=[]
		c=0
		for read in samfile.fetch(chrom, scan_l_pos,scan_r_pos):
			if (read.cigarstring != None) and read_can_shown(read,scan_l_pos,scan_r_pos) and read.mapping_quality >10:
				new_bases,new_base_quality=rearrange_string(read)
				if read_corner_shown(read,scan_l_pos,scan_r_pos,new_bases):
					read_package.append((new_bases,new_base_quality,read))
					c=c+1
		dic=read_to_dictionary(read_package, scan_r_pos,height)
		print ('right: ok before drawn')
		# draw_pgn('right',dic,width,height,scan_l_pos,scan_r_pos,img_name)

		if not c==0 : #图片上有read
			draw_pgn('right',dic,width,height,scan_l_pos,scan_r_pos,img_name)
			imf=open("/mnt/hde/gao/wj/keras/right/homo_right.txt",'a')
			imf.write(tmp[0]+'_'+tmp[1]+'.'+label+'.right.png'+"\t"+label+"\n")
			imf.close()
		# else : # 是一张空白图片，就不画这个图了！
		# 	imf=open("/mnt/hde/gao/wj/SV_image/r_blank_region.txt",'a')
		# 	imf.write(tmp[0]+"\t"+tmp[1]+"\t"+tmp[2]+"\t"+label+"\n")
		# 	imf.close()

		# samfile.close()

		# c=c+1
		# if c==3:
		# 	break
		# break


if __name__ == '__main__':
    main()
