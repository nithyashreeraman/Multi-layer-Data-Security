# -*- coding: utf-8 -*-
import time
time.clock = time.time
import base64
with open("images/processed.png", "rb") as img_file:
    BI = base64.b64encode(img_file.read())

BI = BI.decode("utf-8")

import hashlib
import cv2
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from AESVisualCryptEncrypt import text,K

SK = hashlib.sha256(K.encode())

print("The hexadecimal equivalent of SHA256 is : ")
print(SK.hexdigest())
"""# Decryption"""

cipher = text

P = cv2.imread('shares/P.png')
R = cv2.imread('shares/R.png')

h = np.shape(P)[0]
w = np.shape(P)[1]

CK = np.ones((h,w,1), dtype = 'uint8')

for i in range(h):
    for j in range(w):
        ck = P[i][j][0] ^ R[i][j][0]
        CK[i][j][0] = ck

K1 = []
for i in range(len(CK)):
    K1.append(0)

for i in range(len(CK)):
    count = 0
    for j in range(len(CK[i])):
        if CK[i][j][0] == 0:
            count += 1
    K1[i] = chr(count)

K1 = "".join(K1)

K1

SK1 = hashlib.sha256(K1.encode())

print("The hexadecimal equivalent of SHA256 is : ")
print(SK1.hexdigest())

from Crypto.Cipher import AES
from hashlib import sha256
from base64 import b64decode

class AESCipher:
    def __init__(self,data,key):
        self.block_size = 16
        self.data = data
        self.key = sha256(key.encode()).digest()[:32]
        self.pad = lambda s: s + (self.block_size - len(s) % self.block_size) * chr (self.block_size - len(s) % self.block_size)
        self.unpad = lambda s: s[:-ord(s[len(s) - 1:])]

    def decrypt(self):
        cipher_text = b64decode(self.data.encode())
        iv = cipher_text[:self.block_size]
        cipher = AES.new(self.key,AES.MODE_OFB,iv)
        return self.unpad(cipher.decrypt(cipher_text[self.block_size:])).decode()

xdf = pd.DataFrame(columns = ['1','2'])
a = []
b = []
for i in P:
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

LRmodel = LinearRegression()
LRmodel.fit(xdf,ydf)

zdf = pd.DataFrame(columns = ['1','2'])
a = []
b = []
for i in CK:
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

cipher = cipher.split(' ')

txt = []
for each in cipher:
    try:
        ch = ord(each) - x + y
        txt.append(int(ch))
    except:
        print(each)

text = ""
for t in txt:
    text += chr(t)

de = AESCipher(text,SK1.hexdigest()).decrypt()


de = de.encode("utf-8")

with open("shares/DecryptedImg.png", "wb") as fh:
    fh.write(base64.decodebytes(de))

print("Image is saved 'DecryptedImg.png' ...")
