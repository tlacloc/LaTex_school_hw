#!/usr/bin/env python
# coding: utf-8

# # Instrucciones
# El objetivo de este ejercicio es resolver de manera numérica la siguiente integral
# 
# \\[\vec{F}_{C'\rightarrow C}= \frac{\mu_0}{4\pi}\int_0^{2\pi}\int_0^{2\pi} a \left(-\sin\varphi\hat{i} + a\cos\varphi\hat{j}\right)\text{d}\varphi' \times \left[ a \left(-\sin\varphi'\hat{i} + a\cos\varphi'\hat{j}\right)\text{d}\varphi'\times \frac{(\vec{r} - \vec{r}')}{|\vec{r} - \vec{r}'|^3}\right]\\]
# 
# Donde $\vec{r} - \vec{r'}$ esta dado por
# 
# \\[\vec{r} - \vec{r'}=a(\cos\varphi-\cos\varphi')\hat{i} + a(\sin\varphi-\sin\varphi')\hat{j} + d\hat{k}\\]
# 
# Y entonces $|\vec{r} - \vec{r}'|$ queda de la siguiente forma
# 
# \\[|\vec{r} - \vec{r}'| = \sqrt{a^2(\cos\varphi-\cos\varphi')^2 + a^2(\sin\varphi-\sin\varphi')^2 + d^2}\\]
# 
# Dicha integral se trabajo con el fin de hacerla más fácil expresar en el programa, de modo que se separó en tres integrales, las cuales son:
# 
# \\[A_1 = \frac{\mu_0a^2}{4\pi}\int_0^{2\pi}\int_0^{2\pi}\frac{ - a\cos\varphi\left[\cos{(\varphi-\varphi')}-1\right]\hat{i} }{\left(a^2(\cos\varphi-\cos\varphi')^2 + a^2(\sin\varphi-\sin\varphi')^2 + d^2\right)^{3/2}}\text{d}{\varphi}\text{d}{\varphi'}\\]
# 
# \\[A_2 = \frac{\mu_0a^2}{4\pi}\int_0^{2\pi}\int_0^{2\pi}\frac{ - a\sin\varphi\left[\cos{(\varphi-\varphi')}-1\right]\hat{j} }{\left(a^2(\cos\varphi-\cos\varphi')^2 + a^2(\sin\varphi-\sin\varphi')^2 + d^2\right)^{3/2}}\text{d}{\varphi}\text{d}{\varphi'}\\]
# 
# \\[A_3 = \frac{\mu_0a^2}{4\pi}\int_0^{2\pi}\int_0^{2\pi}\frac{ - d\cos{(\varphi - \varphi')}\hat{k}  }{\left(a^2(\cos\varphi-\cos\varphi')^2 + a^2(\sin\varphi-\sin\varphi')^2 + d^2\right)^{3/2}}\text{d}{\varphi}\text{d}{\varphi'}\\]

# # Cargando librerías

# In[1]:


import sys
sys.path.append('numeric_integration')
from numeric_integration.integration import Integrate
from numpy import sin, cos, pi


# # Definiendo la función

# In[2]:


def test(x, y):
    return sin(x)*sin(y)


# In[3]:


# definiendo el valor de las constantes de la integral:
mu_0 = 1.25663706212*(10**(-6))
a = 1
d = 0.01
pres = 500 #mientras mas grande, mejora la precisión del resultado (aumenta le tiempo de calculo)


# In[4]:


def A1(x, y):
    return -(a*cos(x)*(cos(x-y)-1))/(((a**2)*(cos(x)-cos(y))**2 + (a**2)*(sin(x)-sin(y))**2 + d**2)**(3/2))


# In[5]:


def A2(x, y):
    return -(a * sin(x)*(cos(x-y)-1))/(((a**2)*(cos(x)-cos(y))**2 + (a**2)*(sin(x)-sin(y))**2 + d**2)**(3/2))


# In[6]:


def A3(x, y):
    return -(d*cos(x-y))/(((a**2)*(cos(x)-cos(y))**2 + (a**2)*(sin(x)-sin(y))**2 + d**2)**(3/2))


# # Resultado

# ## A1

# In[7]:


# Build an Integrate object
integral = Integrate(A1)

# Calculate the integral
A1_result = integral.double_integral([[0, 2*pi], [0, 2*pi]], precision=pres)

# Show the result
print("The result is", A1_result)


# ## A2

# In[8]:


# Build an Integrate object
integral = Integrate(A2)

# Calculate the integral
A2_result = integral.double_integral([[0, 2*pi], [0, 2*pi]], precision=pres)

# Show the result
print("The result is", A2_result)


# ## A3

# In[9]:


# Build an Integrate object
integral = Integrate(A3)

# Calculate the integral
A3_result = integral.double_integral([[0, 2*pi], [0, 2*pi]], precision=pres)

# Show the result
print("The result is", A3_result)


# ## Suma de las integrales

# In[10]:


#componentes

iv = ((mu_0*a**2)/(4*pi))*A1_result
jv = ((mu_0*a**2)/(4*pi))*A2_result
kv = ((mu_0*a**2)/(4*pi))*A3_result


# ## Fuerza total

# In[11]:


print('F =', iv , 'i +', jv , 'j +', kv , 'k'  )


# ## Modulo de la fuerza

# In[12]:


print('|F| =', iv**2 + jv**2 + kv**2)


# In[ ]:




