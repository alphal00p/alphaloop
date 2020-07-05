
#include <complex.h>
#include <signal.h>

int evaluate_0(double complex[], int conf, double complex*);
int get_rank_0(int conf);

int evaluate(double complex lm[], int diag, int conf, double complex* out) {
    switch(diag) {
		case 0: return evaluate_0(lm, conf, out);
		default: raise(SIGABRT);
    }
}

int get_buffer_size() {
    return 2;
}

int get_rank(int diag, int conf) {
    switch(diag) {
		case 0: return get_rank_0(conf);
		default: raise(SIGABRT);
    }
}
