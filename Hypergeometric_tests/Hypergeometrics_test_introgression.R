#take numbers from Venn diagram

#phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
#where q=size of overlap-1; m=number of upregulated genes in experiment #1; 
#n=(total number of genes on platform-m); k=number of upregulated genes in experiment #2.


#‘pick m balls from a jar of n balls, then repeat by picking k balls: what is the significance of exactly q overlap. 

###################
#Mias#
###################
#introgression and selection
q=98+19+4+10-1
m=401+30+1+1+98+19+4+10
k=561+98+65+19+4+4+7+10
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#8.672533e-89

#introgression and transcriptomics
q=30+19+4+1-1
m=401+30+1+1+98+19+4+10
k=3052+64+65+4+19+4+30+1
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.7666546

#introgression and CNV
q=4+1+1+10-1
m=401+30+1+1+98+19+4+10
k=236+64+4+7+4+10+1+1
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.0003580077

#selection and transcriptomics
q=65+19+4+4-1
m=561+98+65+19+4+4+7+10
k=3052+64+65+4+19+4+30+1
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.08749724

#selection and CNV
q=4+4+10+7-1
m=561+98+65+19+4+4+7+10
k=236+64+4+7+4+10+1+1
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#7.476793e-07

#transcriptomics and CNV
q=64+4+4+1-1
m=3052+64+65+4+19+4+30+1
k=236+64+4+7+4+10+1+1
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#2.558926e-10


###################
#Selection#
###################
#Mias and Buko
q=36+16+133+69-1
m=252+36+80+16+95+133+87+69
k=280+36+16+22+133+19+69+100
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#1.712507e-239

#Mias and Piek
q=80+16+95+133-1
m=252+36+80+16+95+133+87+69
k=278+48+80+95+16+133+22+19
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0

#Mias and Kato
q=95+133+87+69-1
m=252+36+80+16+95+133+87+69
k=308+48+95+87+133+69+19+100
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0

#Buko and Piek
q=16+22+133+19-1
m=280+36+16+22+133+19+69+100
k=278+48+80+95+16+133+22+19
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#1.064784e-157

#Kato and Piek
q=48+95+133+19-1
m=308+48+95+87+133+69+19+100
k=278+48+80+95+16+133+22+19
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#1.007099e-284

#Buko and Kato
q=133+19+100+69-1
m=280+36+16+22+133+19+69+100
k=308+48+95+87+133+69+19+100
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0

###################
#Introgression#
###################
#Mias and Buko
q=36+24+63+40-1
m=228+36+48+24+62+63+63+40
k=630+36+24+69+63+27+40+78
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#3.252753e-111

#Mias and Piek
q=48+24+62+63-1
m=228+36+48+24+62+63+63+40
k=531+64+48+62+24+63+69+27
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#1.990556e-161

#Mias and Kato
q=62+63+63+40-1
m=228+36+48+24+62+63+63+40
k=445+64+62+63+63+40+27+78
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#1.543964e-211

#Buko and Piek
q=24+69+63+27-1
m=630+36+24+69+63+27+40+78
k=531+64+48+62+24+63+69+27
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#1.796117e-97

#Kato and Piek
q=64+62+63+27-1
m=445+64+62+63+63+40+27+78
k=531+64+48+62+24+63+69+27
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#1.315466e-145

#Buko and Kato
q=63+27+40+78-1
m=630+36+24+69+63+27+40+78
k=445+64+62+63+63+40+27+78
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#2.841984e-128


###################
#CNV#
###################
#Mias and Buko
q=12+16+70+21-1
m=96+12+37+16+42+70+33+21
k=60+12+16+8+70+10+21+22
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#1.077945e-182

#Mias and Piek
q=37+16+42+70-1
m=96+12+37+16+42+70+33+21
k=117+32+37+42+16+70+8+10
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#9.089772e-251

#Mias and Kato
q=42+70+33+21-1
m=96+12+37+16+42+70+33+21
k=108+32+42+33+70+21+10+22
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#2.823872e-251

#Buko and Piek
q=16+8+70+10-1
m=60+12+16+8+70+10+21+22
k=117+32+37+42+16+70+8+10
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#2.308261e-149

#Kato and Piek
q=32+42+70+10-1
m=108+32+42+33+70+21+10+22
k=117+32+37+42+16+70+8+10
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#2.428716e-223

#Buko and Kato
q=70+10+21+22-1
m=60+12+16+8+70+10+21+22
k=108+32+42+33+70+21+10+22
n=31072-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#1.478541e-189







