# HomomorphicFiltering
Homomorphic filtering using opencv

This code is based on paper Butterworth equations for homomorphic ltering of images Computers in Biology and Medicine 28 (1998) 169±181

The modified Butterworth flter function is implemented  as follows for an image with rows line and cols column :

r=((i/rows)^2+(j/cols)^2)^(1/2) where i,j are frequency

Butterworh = (1-1/(1+(r/a)^n))*d+e

Value of r, a, d, n, and e can be modified using upper or lower case r a d n e
