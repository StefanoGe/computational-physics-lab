// test

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

typedef struct {
    float *val;
    int length;
} ArrayFloat;

// Funzione per allocare un array di float
ArrayFloat allocArrF(int length) {
    ArrayFloat array;
    array.length = length;
    array.val = (float *)malloc(length * sizeof(float));  // Assicurati di includere stdlib.h per malloc
    if (array.val == NULL) {
        printf("Errore nell'allocazione della memoria\n");
        exit(1);  // Uscita in caso di errore
    }
    return array;
}

ArrayFloat buildArrF(int length, ...) {
    va_list list;
    va_start(list, length);

    // Allocazione dell'array di float
    ArrayFloat array = allocArrF(length);

    for (int i = 0; i < length; i++) {
        // L'argomento passato a va_arg è trattato come double, quindi facciamo il cast a float
        float temp = (float)va_arg(list, double);
        array.val[i] = temp;
        printf("Valore nell'array: %f\n", array.val[i]);
        printf("Valore in va_arg: %f\n", temp);
    }

    va_end(list);
    return array;
}

int main() {
    // Esempio di utilizzo della funzione variadica
    ArrayFloat array = buildArrF(3, 3.14f, 2.71f, 1.61f);

    // Stampa dei valori nell'array
    for (int i = 0; i < array.length; i++) {
        printf("Array[%d]: %f\n", i, array.val[i]);
    }

    // Non dimenticare di liberare la memoria allocata!
    free(array.val);

    return 0;
}
