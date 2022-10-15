setwd("/Users/a.i/StatisticalProgramming")
a <- scan("pg10.txt",what="character",skip=104)
n <- length(a)
a <- a[-((n-2886):n)]
a <- a[-grep("[0123456789]:[0123456789]",a)]

#(a)
unique_words = unique(tolower(a))
unique_words

#(b)
index_vector = match(tolower(a),unique_words)
index_vector

#(c)
tab = tabulate(index_vector)
tab

#(d)
head(sort(tab,decreasing=TRUE),n=500)
###the word that occurs 144 times is the least of the 500 most common words.

#(e)
b = (unique_words[tab>=144])
b
