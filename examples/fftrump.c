#include "hawaii.h"

double pes(double *q) { UNUSED(q); return 0; }
void dpes(double *q, double *dVdq) { UNUSED(q); UNUSED(dVdq); }

int main()
{
    String_Builder sb = {0};
    CFnc cf = {0};

    if (!read_correlation_function("examples/CF-F-300.0.txt", &sb, &cf)) {
        printf("ERROR: could not read the file!\n");
    } 
   
    printf("%s", sb.items);

    printf("len = %zu\n", cf.len);
    printf("t[0] = %lf\n", cf.t[0]);
    printf("t[-1] = %lf\n", cf.t[cf.len - 1]);

    return 0;  
}
