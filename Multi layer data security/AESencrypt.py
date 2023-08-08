from AESencryptfunc import * #import AESencryptfunc module to use functions created for this program
import math #import math module to use function such as ceiling

_PLAIN_TEXT_PATH = "plainText.txt"
_CIPHER_TEXT_PATH = "Input.txt"

# Each round comprises of 4 steps :
#
# SubBytes
# ShiftRows
# MixColumns
# Add Round Key
# The last round doesn’t have the MixColumns round.
# set passphrase to be a 16 characters, 16 characters * 8 bits = 128 bits strength

PassPhrase=""
while(len(PassPhrase)!=16):
    print("Enter in the 16 character passphrase to encrypt your text file %s" %_PLAIN_TEXT_PATH)
    PassPhrase=input()#takes in user input of char, eg. "Iwanttolearnkung"
    if(len(PassPhrase)<16):#check if less than 16 characters, if so add one space character until 16 chars
        while(len(PassPhrase)!=16):
            PassPhrase=PassPhrase+"\00"
    if(len(PassPhrase)>16):#check if bigger than 16 characters, if so then truncate it to be only 16 chars from [0:16]
        print("Your passphrase was larger than 16, truncating passphrase.")
        PassPhrase=PassPhrase[0:16]

#open plaintext.txt file to read and encrypt
file=open(_PLAIN_TEXT_PATH, "r")
message=(file.read())
print("Inside your plaintext message is:\n%s\n" % message)
file.close()

message=BitVector(textstring=message)
message=message.get_bitvector_in_hex()
replacementptr=0
while(replacementptr<len(message)):
    if(message[replacementptr:replacementptr+2]=='0a'):
        message=message[0:replacementptr]+'0d'+message[replacementptr:len(message)]
        replacementptr=replacementptr+4
    else:
        replacementptr=replacementptr+2

message=BitVector(hexstring=message)
message=message.get_bitvector_in_ascii()
#set up some parameters
start=0#set starting pointer for the part to encrypt of the plaintext
end=0#set ending pointer for the part to encrypt of the plaintex
length=len(message)#check the entire size of the message
loopmsg=0.00#create a decimal value
loopmsg=math.ceil(length/16)+1#use formula to figure how long the message is and how many 16 character segmentss must be encrypted
outputhex=""#setup output message in hex

#A Key Schedule algorithm is used to calculate all the round keys from the key.
# #So the initial key is used to create many different round keys which will be used
#in the corresponding round of the encryption.
#need to setup roundkeys here
PassPhrase=BitVector(textstring=PassPhrase)
roundkey1=findroundkey(PassPhrase.get_bitvector_in_hex(),1)
roundkey2=findroundkey(roundkey1,2)
roundkey3=findroundkey(roundkey2,3)
roundkey4=findroundkey(roundkey3,4)
roundkey5=findroundkey(roundkey4,5)
roundkey6=findroundkey(roundkey5,6)
roundkey7=findroundkey(roundkey6,7)
roundkey8=findroundkey(roundkey7,8)
roundkey9=findroundkey(roundkey8,9)
roundkey10=findroundkey(roundkey9,10)
roundkeys=[roundkey1,roundkey2,roundkey3,roundkey4,roundkey5,roundkey6,roundkey7,roundkey8,roundkey9,roundkey10]

#set up FILEOUT to write
FILEOUT = open(_CIPHER_TEXT_PATH, 'w')

# set up the segement message loop parameters
for y in range(1, loopmsg): # loop to encrypt all segments of the message
    if(end+16<length): #if the end pointer is less than the size of the message, then set the segment to be 16 characters
        plaintextseg = message[start:end + 16]
    else: #or else if the end pointer is equal to or greator than the size of the message
        plaintextseg = message[start:length]
        for z in range(0,((end+16)-length),1): #run a while loop to pad the message segement to become 16 characters, if it is 16 already the loop will not run
            plaintextseg = plaintextseg+"\00"
            #plaintextseg2=BitVector(textstring=plaintextseg)
            #print(plaintextseg2.get_bitvector_in_hex())

    #add round key zero/ find round key one
    bv1 = BitVector(textstring=plaintextseg)
    bv2 = PassPhrase
    resultbv=bv1^bv2
    myhexstring = resultbv.get_bitvector_in_hex()


    for x in range(1, 10):  # loop through 9 rounds
        # sub byte
        #In this step each byte is substituted by another byte. Its performed using a lookup table also called the S-box.
        # This substitution is done in a way that a byte is never substituted by itself
        # and also not substituted by another byte which is a compliment of the current byte.
        # The result of this step is a 16 byte (4 x 4 ) matrix like before.
        myhexstring = resultbv.get_bitvector_in_hex()
        temp1=subbyte(myhexstring)

        # shift rows
        #Each row is shifted a particular number of times.
#The first row is not shifted
#The second row is shifted once to the left.
#The third row is shifted twice to the left.
        temp2=shiftrow(temp1)

        # mix column
        #This step is basically a matrix multiplication.
        # Each column is multiplied with a specific matrix
        # and thus the position of each byte in the column is changed as a result.
        #This step is skipped in the last round.
        bv3 = BitVector(hexstring=temp2)
        newbvashex=mixcolumn(bv3)
        newbv=BitVector(hexstring=newbvashex)

        #add roundkey for current round
        #Now the resultant output of the previous stage is XOR-ed with the corresponding round key.
        #Here, the 16 bytes is not considered as a grid but just as 128 bits of data.
        #After all these rounds 128 bits of encrypted data is given back as output.
        #This process is repeated until all the data to be encrypted undergoes this process.
        bv1 = BitVector(bitlist=newbv)
        bv2 = BitVector(hexstring=roundkeys[x-1])
        resultbv = bv1 ^ bv2
        myhexresult = resultbv.get_bitvector_in_hex()

    #start round 10
    # sub byte round 10
    myhexstring = resultbv.get_bitvector_in_hex()
    temp1=subbyte(myhexstring)

    # shift rows round 10
    temp2=shiftrow(temp1)

    # add round key round 10
    newbv = BitVector(hexstring=temp2)
    bv1 = BitVector(bitlist=newbv)
    bv2 = BitVector(hexstring=roundkeys[9])
    resultbv = bv1 ^ bv2
    myhexstring = resultbv.get_bitvector_in_hex()

    #set encrypted hex segement of message to output string
    outputhextemp = resultbv.get_hex_string_from_bitvector()
    FILEOUT.write(outputhextemp)
    start = start + 16 #increment start pointer
    end = end + 16 #increment end pointer

# encrypted output hex string to specified cipherhex file
FILEOUT.close()

file2=open(_CIPHER_TEXT_PATH, "r")
print("The output hex value for the entire message is:\n%s\n" % file2.read())
file2.close()