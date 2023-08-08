# -*- coding: utf-8 -*-
import time
time.clock = time.time
import base64
with open("images/processed.png", "rb") as img_file:
    BI = base64.b64encode(img_file.read())

BI = BI.decode("utf-8")

import hashlib
print("Enter in a key for visual encryption");
K=input()

SK = hashlib.sha256(K.encode())

print("The hexadecimal equivalent of SHA256 is : ")
print(SK.hexdigest())



# AES 256 in OFB mode:
from Crypto.Cipher import AES
from Crypto.Random import new as Random
from hashlib import sha256
from base64 import b64encode,b64decode

class AESCipher:
    def __init__(self,data,key):
        self.block_size = 16
        self.data = data
        self.key = sha256(key.encode()).digest()[:32]
        self.pad = lambda s: s + (self.block_size - len(s) % self.block_size) * chr (self.block_size - len(s) % self.block_size)
        self.unpad = lambda s: s[:-ord(s[len(s) - 1:])]

    def encrypt(self):
        plain_text = self.pad(self.data)
        iv = Random().read(AES.block_size)
        cipher = AES.new(self.key,AES.MODE_OFB,iv)
        return b64encode(iv + cipher.encrypt(plain_text.encode())).decode()

# Encrypting image using base 64 encoded text and hashed key - SHA256
# AES-256
c = AESCipher(BI,SK.hexdigest()).encrypt()

"""Additional Encryption"""

import numpy as np
import cv2
w = 255
h = len(K)
# creating new Image C of size(w,h)
# initializing as blank
C = np.ones((h,w,1), dtype = 'uint8')

# Filling pixels in C
for i in range(h):
    j = ord(K[i])
    for k in range(w):
        if k < j:
            C[i][k][0] = 0
        else:
            break

# print(C[0])

# Dividing C into R and P
# initializing R and P of same size as C
R = np.ones((h,w,3), dtype = 'uint8')
P = np.ones((h,w,3), dtype = 'uint8')

for i in range(h):
    for j in range(w):
        r = np.random.normal(0,1,1)
        #print("r: ", r)
        R[i][j][0] = r

for i in range(h):
    for j in range(w):
        p = R[i][j][0] ^ C[i][j][0]
        P[i][j][0] = p

filename = 'shares/R.png'
cv2.imwrite(filename, R)

filename = 'shares/P.png'
cv2.imwrite(filename, P)

"""Encrypt CipherText further"""

import pandas as pd

xdf = pd.DataFrame(columns = ['1','2'])
a = []
b = []
cnt = 0
for i in P:
    #print("I: ", i)
    #print("I Shape: ", i.shape)
    cnt+=1
    k = 0
    n1 = []
    n2 = []
    for j in i:
        #print("J SHape: ", j.shape)
        if k%2==0:
            n1.append(np.sum(j))
        else:
            n2.append(np.sum(j))
        k += 1
    a.append(sum(n1))
    b.append(sum(n2))
# print("Cnt: ", cnt)
xdf['1'] = a
xdf['2'] = b

ydf = pd.DataFrame(columns = ['1','2'])
a = []
b = []
for i in R:
    k = 0
    n1 = []
    n2 = []
    for j in i:
        if k%2==0:
            n1.append(np.sum(j))
        else:
            n2.append(np.sum(j))
        k += 1
    a.append(sum(n1))
    b.append(sum(n2))
ydf['1'] = a
ydf['2'] = b

print(xdf)

print(ydf)

"""Fit LR Model

"""

from sklearn.linear_model import LinearRegression
LRmodel = LinearRegression()
LRmodel.fit(xdf,ydf)

zdf = pd.DataFrame(columns = ['1','2'])
a = []
b = []
for i in C:
    k = 0
    n1 = []
    n2 = []
    for j in i:
        if k%2==0:
            n1.append(np.sum(j))
        else:
            n2.append(np.sum(j))
        k += 1
    a.append(sum(n1))
    b.append(sum(n2))
zdf['1'] = a
zdf['2'] = b


predict = LRmodel.predict([[sum(zdf['1']),sum(zdf['2'])]])


x = round(predict[0][0])%26
y = round(predict[0][1])%26

txt = []
for each in c:
    ch = ord(each) + x - y
    txt.append(int(ch))

len(txt)

text = ""
for t in txt:
    text += chr(t) + " "

text

f = open("shares/cipher.txt",'a',encoding='utf-8')
f.write(text)
f.close()