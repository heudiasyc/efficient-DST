import os

f = open('./license-header','r')
header = f.read()
f.close()

for file_path in os.listdir("./"):
	full_path = os.path.join('.', file_path)
	if os.path.isfile(full_path):
		print("prepending file ", full_path)
		f = open(full_path,'r')
		temp = f.read()
		f.close()

		f = open(full_path, 'w')
		f.write(header)
		f.write('\n')
		f.write(temp)
		f.close()
