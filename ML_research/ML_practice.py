# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 23:19:41 2018

@author: aalco
"""


from sklearn import svm, datasets
import matplotlib.pyplot as plt
digits = datasets.load_digits()

clf = svm.SVC(gamma = 0.01, C = 100)

x,y = digits.data[:-1] , digits.data[:-1]

clf.fit(x,y)


print(clf.predict(digits.data[-1]))
plt.imshow(digits.images[-1])