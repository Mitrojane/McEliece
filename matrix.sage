load('ikk_simulation.py')
n = 5; R = GF(2);
for i in range(2, 15):
   n = 2*i;
   v = random_vector(R, n);
   print('size:')
   print(n);

   t0 = time.time()
   m = random_nonsingular_matrix(R, n, n);
   t1 = time.time()
   print("gen_random_matrix: %s s" % (t1-t0))

   t0 = time.time()
   m2=~(m)
   t1 = time.time()
   print("inverse_matrix: %s s" % (t1-t0))

   t0 = time.time()
   tmp = m*m2;
   t1 = time.time()
   print("multiply_matrix: %s s" % (t1-t0))




