from PIL import Image
import math
import time


start_time = time.time()
_TEXT_PATH = "Input.txt"
_IMAGE_PATH = "images/original.png"
_SAVE_PATH = "images/processed.png"
_PERFECT_RATIO = True


#Handling the message to be encrypted
#1. Reading the text inside file specified in _TEXT_PATH
#2. Convert the message into binary ascii splitted by 4
#   For example 01011101 will be appended to an array as [01,01,11,01]
#   And that goes for each letter 
#3. Return the array
def readFile(path):

    fileR = open(path, 'r')
    msg = fileR.read()

    binaArr = []
    for element in msg:  
        binaArr.append(format(ord(element), '08b'))
    splitBinaArray = []
    for x in binaArr:
        splitBinaArray.append(x[0:2])
        splitBinaArray.append(x[2:4])
        splitBinaArray.append(x[4:6])
        splitBinaArray.append(x[6:8])
    
    for x in range(4):
        if (len(splitBinaArray)%3!=0):
            splitBinaArray.append('00')
        else:
            break

    return splitBinaArray


#Getting text as binary
bina = readFile(_TEXT_PATH)



#Opening the original image to be processed and making a copy of it to save later.
originalImage = Image.open(_IMAGE_PATH)
processedImage = originalImage.copy()

#Prining some information
imgMode = processedImage.mode
pixelNeeded = int(len(bina)/3)
pixelAvailable = processedImage.size[0]*originalImage.size[1]
print("Image mode", imgMode)
print("Pixles Needed = ",pixelNeeded)
print("Pixles Avilable",pixelAvailable)


#Stop the application if the picture is smaller than text
#   1 CHAR = 4 RGB
#   3 RGB  = 1 PIXEL
#   NEW LINE = 2 CHAR = 8 RGB
if(pixelNeeded>pixelAvailable):
   
    new_x =0
    new_y =0
    x=(pixelNeeded/pixelAvailable)

    if(_PERFECT_RATIO):
        r=math.ceil(math.sqrt(x))
        new_x =processedImage.size[0]*r
        new_y =originalImage.size[1]*r
    else:
        r=math.sqrt(x)
        new_x =math.ceil(processedImage.size[0]*r)
        new_y =math.ceil(originalImage.size[1]*r)
        
    
    
    print("!!! IMAGE HAS BEEN SCALED UP TO MAINTAIN EXTRA DATA !!!")
    print("!!! OLD SIZE",processedImage.size[0],"x",originalImage.size[1],"=",processedImage.size[0]*processedImage.size[1],"  RATIO : ",processedImage.size[0]/processedImage.size[1]," !!!")
    print("!!! NEW SIZE",new_x,"x",new_y,"=",new_x*new_y,"   RATIO : ",new_x/new_y," !!!")
    processedImage=processedImage.resize((new_x,new_y),Image.ANTIALIAS)
    

    



#Looping through the pixels and replacing the last 2 bits of each color channel
#with an element of the binary array for example Array = [11,10,01,10,01]
#  (RGB1)=(0000000,00000000,0000000) (RGB2)=(0000000,00000000,0000000)
#  After proccessing
#  (RGB1)=(0000011,00000010,0000001) (RGB2)=(0000010,00000001,0000000)
i=0
halt=False
for x in range(processedImage.size[0]):
    for y in range(processedImage.size[1]):
        if(i<len(bina)):  
            _rgb = processedImage.getpixel((x,y))
            _r=int(format(_rgb[0],'08b')[0:6]+bina[i],2)
            _g=int(format(_rgb[1],'08b')[0:6]+bina[i+1],2)
            _b=int(format(_rgb[2],'08b')[0:6]+bina[i+2],2)
            
            processedImage.putpixel((x,y),(_r,_g,_b))
            #print(processedImage.getpixel((x,y)))
            i=i+3
        else:
            halt=True
            break
    if(halt):
        break




#Saving the processed image with full quality to avoid loss of data.
processedImage.save(_SAVE_PATH,format="png",quality=100)
print("ENCODED IN IMAGE => ",_SAVE_PATH)
print("ENCODED TOOK %s seconds" % (time.time() - start_time))


