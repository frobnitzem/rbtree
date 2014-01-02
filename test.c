#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "rbtree.h"

struct dirent;
struct dirent {
    int n;
    unsigned char mark;
    struct dirent *L, *R;
} dex; // member for addressing purposes only

// could also have used NULL
static struct dirent nil = {-1, 0, NULL, NULL};

static int int_cmp(const void *ai, const void *bi) {
    const struct dirent *a = ai;
    const struct dirent *b = bi;
    return a->n - b->n;
}

static rbop_t rbinf = {
    .cmp = int_cmp,
    .coff = (void *)&(dex.L) - (void *)&dex,
    .boff = (void *)&(dex.mark) - (void *)&dex,
    .nil = &nil,
    .mask = 1, // use smallest bit
};

//static int N = 16;
static int N = 1 << (4*3); // (~ 4k)
//static int N = 1 << (4*3+10); // (~ 4k)*1024

static void dot_rec(FILE *f, struct dirent *a, int n);
void tree_to_dot(FILE *f, struct dirent *a);
int show_tree(char *name, struct dirent *a, int waitfor);

int main(int argc, char **argv) {
    int i, j, k;
    int ord[N];
    char buf[32];
    struct dirent *ent, *ret;
    //struct dirent *tree = NULL;
    void *tree = rbinf.nil; // actually struct dirent *
    // will segfault is tree is NULL (unless rbinf.nil == NULL)
    size_t len;
    
    if( (ent = malloc(N*sizeof(struct dirent))) == NULL) {
        perror("malloc");
        return 2;
    }
    if(argc >= 2) {
        i = atoi(argv[1]);
        srand(i);
        printf("Seeded rand with %d\n", i);
    }
    for(i=0; i<N; i++) {
        ent[i].n = i;
        ord[i] = i;
    }
    printf("Testing %d additions.\n", N);
    for(j=0; j<N; j++) {
        k = random() % (N-j) + j; // randomize ord in-place
        i = ord[k];
        ord[k] = ord[j];
        ord[j] = i;
        //printf("Adding %d.\n", i);
        if(add_node(&tree, ent+i, &rbinf) != rbinf.nil) goto err;
        /*if(j < 10) {
            show_tree("test.dot", tree, 1);
        }*/
    }
    printf("Finished addition phase.\n");
    //show_tree("test.dot", tree, 0);

    printf("Testing false del.\n");
    i = -1; // non-existent node
    if( (ret = del_node(&tree, (void *)&i, &rbinf)) != rbinf.nil) {
        printf("Got: %d\n", ret->n);
        goto err;
    }

    printf("Testing %d deletions.\n", N);
    for(j=0; j<N; j++) {
        k = random() % (N-j) + j; // remaining item id.
        i = ord[k];
        ord[k] = ord[j];
        ord[j] = i;
        //printf("\nDeleting %d [OK]: ", i);
        //fgets(buf, sizeof(buf), stdin);
        if( (ret = del_node(&tree, (void *)&i, &rbinf)) == rbinf.nil)
            goto err;
        //printf("Got: %d\n", ret->n);
        /*if(N-j < 10) {
            show_tree("test.dot", tree, 1);
        }*/
    }
    printf("Done!\n");
    printf("Testing false del.\n");
    i = 1;
    if( (ret = del_node(&tree, (void *)&i, &rbinf)) != rbinf.nil) {
        printf("Got: %d\n", ret->n);
        goto err;
    }

    free(ent);
    return 0;

err:
    printf("an error occured.\n");
    free(ent);
    return 1;
}

static void dot_rec(FILE *f, struct dirent *a, int n) {
    fprintf(f, "  %d [", a->n);
    if(get_mask(a, &rbinf)) {
        fprintf(f, "color=\"red\" ");
    } else {
        fprintf(f, "color=\"black\" ");
    }
    fprintf(f, "rank=%d];\n", n);

    if(a->L != &nil) {
        fprintf(f, "%d -> %d;\n", a->n, a->L->n);
        dot_rec(f, a->L, n+1);
    }
    if(a->R != &nil) {
        fprintf(f, "%d -> %d;\n", a->n, a->R->n);
        dot_rec(f, a->R, n+1);
    }
}

void tree_to_dot(FILE *f, struct dirent *a) {
    fprintf(f, "digraph RBTree {\n");
    if(a != &nil) {
        dot_rec(f, a, 0);
    }
    fprintf(f, "}\n");
}

int show_tree(char *name, struct dirent *a, int waitfor) {
    FILE *f;
    char *buf;
    int stat;
    pid_t pid;

    if( (f = fopen(name, "w")) == NULL) {
        perror("Error opening dot output file");
        return 1;
    }
    tree_to_dot(f, a);
    fclose(f);

    if( (pid = fork()) == 0) {
        execlp("dotty", "dotty", name, NULL);
        exit(1);
        /*if(asprintf(&buf, "dot -Tsvg -o %s.svg %s && inkview %s.svg", name, name, name) < 0) {
            exit(1);
        }
        exit(system(buf));*/
    }
    if(waitfor) {
        waitpid(pid, &stat, WUNTRACED);
    }
    return 0;
}


