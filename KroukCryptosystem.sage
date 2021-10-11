class KroukCryptosystem:


    def __init__(self, n, m, g):
        # Construct the Goppa code
        goppa_code = GoppaCode(n, m, g);
        assert goppa_code.generator_matrix().nrows() != 0, "Generator Matrix is empty.";
        k = goppa_code.generator_matrix().nrows();
        
        self._m_GoppaCode = goppa_code;
        self._g = g;
        self._t = g.degree();


        # Set up the random scrambler matrix
        M = matrix(GF(2), n, [random() < 0.5 for _ in range(n ** 2)]);
        while (rank(M) < n):
            M = matrix(GF(2), n, [random() < 0.5 for _ in range(n ** 2)]);

        # Set up the random scrambler matrix
        W = matrix(GF(2), n, [random() < 0.5 for _ in range(n ** 2)]);
        while (rank(W) < n):
            W = matrix(GF(2), n, [random() < 0.5 for _ in range(n ** 2)]);


        # Set up the random scrambler matrix
        U = matrix(GF(2), n, k,[random() < 0.5 for _ in range(n*k)]);
        while (rank(U) < k):
           U = matrix(GF(2), n, k,[random() < 0.5 for _ in range(n*k)]);

        Cn = U*(goppa_code.generator_matrix());

        # Set up the permutation matrix
        rng = list(range(n));
        P = matrix(GF(2), n);
        for i in range(n):
            p = floor(len(rng) * random());
            P[i, rng[p]] = 1;
            rng = rng[:p] + rng[p + 1:];


        # Set up the random diagonal matrix
        zeros = set([]);
        while len(zeros) != (n - self.max_num_errors()):
            zeros.add(floor(n * random()))

        D = matrix(GF(2), n)
        for i in range(n):
            if (not i in zeros):
                 D[i,i] = 1;



        # Remember these values
        self._G = goppa_code.generator_matrix();
        self._M = M;
        self._P = P;
        self._D = D;
        self._W = W;
        self._Cn = Cn;
        self._PublicKeyE = W*D*(Cn +P)*M;
        self._PublicKeyG = self._G*M;

    # This is a help function which will be useful for encryption.
    def _GetRowVectorWeight(self, n):
        weight = 0;
        for i in range(n.ncols()):
            if n[0, i] == 1:
                weight = weight + 1;
        return weight;




    def Encrypt(self, message):
        # Assert that the message is of the right length
        assert (message.ncols() == self.public_key_G().nrows()), "Message is not of the correct length.";

        # Get an error vector, ensuring that there are exactly t errors.
        err_vec = matrix(1, self.goppa_code().generator_matrix().ncols());
        while (self._GetRowVectorWeight(err_vec) < self.max_num_errors()):
            err_vec[0, randint(1, self.goppa_code().generator_matrix().ncols() - 1)] = 1;
        code_word = message * self.public_key_G();
        return copy(code_word + err_vec*self.public_key_E());





    def Decrypt(self, received_word):
        assert (received_word.ncols() == self.public_key_G().ncols()), "Received word is not of the correct length.";

        # Strip off the permutation and decode the received word.
        message = received_word * ~(self._M);
        #print(message);
        tmp = self.goppa_code().Decode(message);
        #print(tmp);
        e_vector = matrix(GF(2), 1,message.ncols());
        for i in range(message.ncols()):
            if(message[0, i] != tmp[0, i]):
                e_vector[0, i] = 1;
        #print(e_vector);
        errors = e_vector * ~(self._P);
        errors = (self._W * self._D).solve_left(errors)

        tmp = received_word - errors * (self.public_key_E());

        # Solve the system to determine the original message.
        message = self.public_key_G().solve_left(tmp);

        return copy(message);



    # Accessors
    def public_key_G(self):
        return copy(self._PublicKeyG);

    # Accessors
    def public_key_E(self):
        return copy(self._PublicKeyE);

    def goppa_code(self):
        return copy(self._m_GoppaCode);

    def max_num_errors(self):
        return copy(self._t);

