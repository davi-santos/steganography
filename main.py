#!/usr/bin/env python
# coding:UTF-8
"""main.py

Usage:
  main.py hide -i <input> -o <output> -f <file> 
  main.py unhide -i <input> -o <output>
  main.py --help

Options:
  -h, --help                Show this help
  -f,--file=<file>          File to hide - arquivo para escondedr
  -i,--in=<input>           Cover Image - imagem de cobertura	
  -o,--out=<output>         Output image (with data embbeded or recovered image) - nome da imagem de saida ou da imagem recuperada
"""
from scipy import misc
from math import ceil
import numpy as np
import docopt

NUM_BITS = 8
MASK = [1,3,7,15]

# --- funções para salvar dados secretos na imagem

def getHeaderSize (kLSB):
	if kLSB == 1: return 39
	if kLSB == 2: return 51
	if kLSB == 3: return 75
	else: return 123

def setCoverImageLSB(coverImageArray1D, kLSB, headerSize, qtyPixels):
	mask = MASK[kLSB-1]
	qtyChannels = qtyPixels * int(ceil(8/kLSB)) + headerSize
	coverImageArray1D[:headerSize] = coverImageArray1D[:headerSize] & (~3)
	coverImageArray1D[headerSize:qtyChannels] = coverImageArray1D[headerSize:qtyChannels] & ~mask

	return coverImageArray1D

def saveDataTo1DVector(secretData, nullVector, kLSB, header):
	
	mask = MASK[kLSB-1]
	step = int(ceil(8/kLSB)) #funciona apenas com k=1, k=2 e k=4

	for item in secretData:
		for i in range(step):
			aux = item & mask
			nullVector[header] = aux
			item = item >> kLSB
			header += 1

	return nullVector

def replaceColorArray(secretDataArray, colorArray, header, qtyChannels):
	i = 0; qtyChannels += header
	for item in secretDataArray[header:qtyChannels]:
		secretDataArray[header] = colorArray[i,secretDataArray[header]]
		i = (i+1) % 3
		header += 1
	return secretDataArray

def joinImages (vetorDadosSecretos, coverImageArray1D):
	return vetorDadosSecretos|coverImageArray1D

def calcKLSB (coverImageMatrix, secretDataImage):
	qtyPixels = coverImageMatrix.size
	qtyData = secretDataImage.size
	canais = [123, 75, 51, 39] 
	count = 1

	for i in range(4):
		qtyCanais = canais.pop()
		sizeCoverImage = qtyPixels - qtyCanais
		if sizeCoverImage == 0: return -1
		kLSB = int(ceil(NUM_BITS*(float(qtyData)/sizeCoverImage)))
		if kLSB == count:
			print("Valor de K na funcao: ", kLSB)
			#if kLSB == 3:
			#	return 4
			return kLSB
		count += 1
	return -1

def createColorArrays (kLSB):
	matrix_size = 2**kLSB
	matrixColors = np.empty((3,matrix_size), dtype ='uint8')
	vector = np.arange(matrix_size, dtype='uint8')

	for i in range(3):
		np.random.shuffle(vector)
		matrixColors[i] = vector
	return matrixColors

def setArrayHeader(kLSB, vectorColors, sizeX, sizeY, emptyArray):
	#mask=3 e k=2 porque usa-se 2 LSBs para o cabeçário
	mask = 3; count = 0; k = 2

	vectorColors = vectorColors.flatten()
	listOfValues = np.concatenate(([kLSB], vectorColors))

	for item in listOfValues:
		for i in range(k):
			emptyArray[count] = item & mask
			item = item >> k
			count += 1

	for item in [sizeX, sizeY]:
		for i in range(12):
			emptyArray[count] = item & mask
			item = item >> k
			count += 1			

	count += 1
	return emptyArray

def setStega (coverImage, secretImage, outImage):
	#open image and create a new one
	coverImage = misc.imread(coverImage)
	secretImage = misc.imread(secretImage)
	sizeX, sizeY, z = secretImage.shape #OK HAHA
	sizeXFinal, sizeYFinal, zFinal = coverImage.shape
	finalImage = np.zeros(coverImage.size, dtype='uint8')

	#change each image to 1D array
	coverImage = coverImage.flatten()
	secretImage = secretImage.flatten()

	#calculate K and header size
	kLSB = calcKLSB(coverImage, secretImage) #OK HAHA
	print("Valor de K atual: ", kLSB)
	#print("valor de k: ", kLSB)
	if kLSB == -1: return 0
	header = getHeaderSize(kLSB) #OK HAHA
	colorArrays = createColorArrays(kLSB) #OK HAHA
	qtyChannels = secretImage.size*int(ceil(8/kLSB))

	#limpar os LSBs da coverImage
	residualImage = setCoverImageLSB(coverImage, kLSB, header,secretImage.size) #OK HAHA ~
	#salvar header na new image
	finalImage = setArrayHeader(kLSB, colorArrays, sizeX, sizeY, finalImage) #OK HAHAHA
	#salvar dados secretos na imagem
	finalImage = saveDataTo1DVector(secretImage, finalImage, kLSB, header)
	#converter dados atraves da matriz
	finalImage = replaceColorArray(finalImage, colorArrays, header, qtyChannels)
	#juntar imagem residual e imagem para esconder
	finalImage = joinImages(finalImage, residualImage)
	finalImage = finalImage.reshape((sizeXFinal, sizeYFinal, zFinal))

	misc.imsave(outImage, finalImage)

# -- funções para recuperar matriz na imagem de cobertura

def recoverHeader(coverImage1DArray):

	#recuperando o valor de K
	mask = 3; shift = 0; valor = 0
	for i in coverImage1DArray[:2]:
		aux = i & mask
		aux = aux << shift
		shift += 2
		valor = valor |  aux

	kLSB = valor
	tamanhoVetores = 3*(2**kLSB)
	vetoresVazio = np.zeros(tamanhoVetores, dtype="uint8")
	canais = 2*3*(2**kLSB) #3 vetores, cada valor ocupa 2 canais, retorna qty canais dos vetores
	canais += 2 #soma dois, pois 2 primeiros canais é K

	shift = 0; valor = 0; cont = 0; index = 0
	for item in coverImage1DArray[2:canais]:
		aux = item & mask
		aux = aux << shift
		shift = (shift+2) % 4
		valor = valor | aux
		cont += 1
		if cont == 2:
			vetoresVazio[index] = valor
			valor = 0
			index += 1
			cont = 0

	limite = canais + 24; cont = 0; valor = 0; index = 0; shift = 0
	dimensao = np.array([0,0])
	
	for item in coverImage1DArray[canais:limite]:
		aux = item & mask
		aux = aux << shift
		shift = (shift+2) % 24
		valor = valor | aux
		cont += 1
		if cont == 12:
			dimensao[index] = valor
			index += 1
			valor = 0
			cont = 0

	vetoresVazio = vetoresVazio.reshape((3,2**kLSB))
	limite += 1
	return kLSB, vetoresVazio, dimensao, limite

def recoverImage (coverImage1DArray, limite, numCanais, kLSB, vectorColors):
	#cria uma nova imagem vazia em 1D
	newImage = np.zeros(numCanais, dtype = "uint8")
	numCanais *= 2
	numCanais += limite
	mask = MASK[kLSB-1]
	aux = 0; i = 0; shift = 0; cont = 0; index = 0; valor = 0

	for item in coverImage1DArray[limite:numCanais]:
		aux = item & mask
		aux, = np.where(vectorColors[i] == aux)[0]
		aux = aux << shift
		valor = valor | aux
		i = (i+1) % 3
		shift += shift+kLSB
		cont += kLSB
		
		if cont == 8:
			newImage[index] = valor
			valor = 0
			index += 1
			cont = 0
			shift = 0
			
	return newImage

def recoverImageAgain (coverImage1DArray, header, numCanais, kLSB, vectorColors, sizeX, sizeY):
	#criar nova Imagem
	tamanho = sizeX*sizeY*3
	newImage = np.zeros((tamanho),dtype="uint8")
	passo = int(ceil(8/kLSB)); mask = MASK[kLSB-1]

	shift = 0; valor = 0; cont = 0

	for oi in range(tamanho):

		for line in range(passo):
			aux = coverImage1DArray[header] & mask
			aux, = np.where(vectorColors[cont] == aux)[0]
			aux = aux << shift
			cont = (cont+1) % 3
			valor = valor | aux
			shift += kLSB
			header += 1

		newImage[oi] = valor
		valor = 0; shift = 0

	return newImage

def getStega (coverImage, outImage):
	#abrir imagem
	coverImage = misc.imread(coverImage)
	coverImage = coverImage.flatten()

	kLSB, vetores, dimensao, limite = recoverHeader(coverImage)
	sizeX, sizeY = dimensao

	channels = int(sizeX*sizeY*3)
	novaImagem = recoverImageAgain (coverImage, limite, channels, kLSB, vetores, sizeX, sizeY)
	novaImagem = novaImagem.reshape((sizeX, sizeY, 3))

	misc.imsave(outImage, novaImagem)

def main():
	args = docopt.docopt(__doc__)
	inputFile = args["--in"]
	outFile = args["--out"]
	
	if args['hide']:
		secretFile = args["--file"]
		setStega (inputFile, secretFile, outFile)
		print("Esteganografia feita com sucesso")

	elif args['unhide']:
		getStega(inputFile, outFile)
		print("Imagem recuperada com sucesso")

if __name__ == "__main__":
	main()
