from PIL import Image
import time


start_time = time.time()
_IMAGE_PATH = "shares/DecryptedImg.png"
_OUTPUT_PATH = "Output.txt"


#Opening the image to be decoded
img = Image.open(_IMAGE_PATH)
print("Image mode", img.mode)
print("Image Size",img.size[0],"x",img.size[1],"=",img.size[0]*img.size[1])


#creating array that holds the hidden bits
bina = []


#Extracting the least siginificant 2bits from each pixel of the image and storing them in the array
for x in range(img.size[0]):
    for y in range(img.size[1]):
        _rgb = img.getpixel((x,y))
        _r_bina = format(_rgb[0],'08b')[6:8]
        _g_bina = format(_rgb[1],'08b')[6:8]
        _b_bina = format(_rgb[2],'08b')[6:8]
        bina.append(_r_bina)
        bina.append(_g_bina)
        bina.append(_b_bina)


#Fixing the length of the array to be devidable by four
#Since ascii is 8 bits and we partitioned it to be 4 parts for each character
#Example character = 11001101 = [11,00,11,01], inside the array
for x in range(4):
    if(len(bina)%4!=0):
        bina.append('00')
    else:
        break

#Creating an array of characters 
characters = []

#We will use this to convert the array into a full string
def convertToString(s):
    new = ""
    for x in s:
        new += x 
    return new
      

#Reading 4 elements in each loop in the binary array
#Then convert them into a character and append it to the array
#
#We stop the array when we have a character that is out of the range of
#ascii letters, we add char=10 because that translate to new line
x=0
while(x<len(bina)):
    char = int(bina[x]+bina[x+1]+bina[x+2]+bina[x+3],2)
    if((char<127 and char>30) or char==10):

        #print(char)
        characters.append(chr(char))
        x=x+4
    else:
        break


text = convertToString(characters[:-1])


#save the extracted string and print in the console

with open(_OUTPUT_PATH, "w") as f:
    f.writelines(str(text))
    #print("=====MESSAGE=====\n",text,"\n=====END MESSAGE=====")
print("SAVED MESSAGE IN => Output.txt")
print("DECODING TOOK %s seconds" % (time.time() - start_time))