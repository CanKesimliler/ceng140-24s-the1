#include <stdio.h>
#define MAX_DNA_LENGTH 3005
#define MAX_PROTEIN_LENGTH 3000

void getInput(char *, char *, int *, int *);
int CheckCode(int key, char (*aminoacids)[6][3], char *readingFrame);
int findExon1(char *, char *, char (*aminoacids)[6][3], int *, int *, int *);
int findExon2(char *, char *, char (*aminoacids)[6][3], int *, int, int , int, int);
int findExon3(char *, char *, char (*aminoacids)[6][3], int *, int *, int *);

int main()
{
    int map[200];
    char aminoacids[21][6][3] = {
        {"CCT", "CCG", "CCA", "CCC", "", ""},       /*gly 0*/
        {"CGT", "CGA", "CGC", "CGG", "", ""},       /*ala 1*/
        {"CAC", "CAT", "CAG", "CAA", "", ""},       /*val 2*/
        {"TAA", "TAG", "TAT", "", "", ""},          /*ile 3*/
        {"AAT", "AAC", "GAG", "GAC", "GAA", "GAT"}, /*leu 4*/
        {"TCA", "TCG", "AGA", "AGG", "AGC", "AGT"}, /*ser 5*/
        {"TGA", "TGG", "TGT", "TGC", "", ""},       /*thr 6*/
        {"GGA", "GGG", "GGT", "GGC", "", ""},       /*pro 7*/
        {"CTA", "CTG", "", "", "", ""},             /*asp 8*/
        {"CTT", "CTC", "", "", "", ""},             /*glu 9*/
        {"TTT", "TTC", "", "", "", ""},             /*lys 10*/
        {"TCT", "TCC", "GCA", "GCG", "GCT", "GCC"}, /*arg 11*/
        {"TTA", "TTG", "", "", "", ""},             /*asn*/
        {"GTT", "GTC", "", "", "", ""},             /*gln*/
        {"ACA", "ACG", "", "", "", ""},             /*cys*/
        {"TAC", "", "", "", "", ""},                /*met*/
        {"ACC", "", "", "", "", ""},                /*trp*/
        {"AAA", "AAG", "", "", "", ""},             /*phe*/
        {"ATA", "ATG", "", "", "", ""},             /*tyr*/
        {"GTA", "GTG", "", "", "", ""},             /*his*/
        {"ATT", "ATC", "ACT", "", "", ""}};         /*stop*/
    int intrn1Start, intrn2Start, intrn1End, intrn2End;
    char dnaSequence[MAX_DNA_LENGTH] = "";
    char proteinSequence[MAX_PROTEIN_LENGTH] = "";
    int proteinLen = -1;
    int dnaLen = -1;
    int ldna, rdna, lpro, rpro;
    intrn1Start = intrn1End = intrn2Start = intrn2End = 0;
    ldna = lpro = 0;
    map[(int)'g'] = 0;
    map[(int)'a'] = 1;
    map[(int)'v'] = 2;
    map[(int)'i'] = 3;
    map[(int)'l'] = 4;
    map[(int)'s'] = 5;
    map[(int)'t'] = 6;
    map[(int)'p'] = 7;
    map[(int)'d'] = 8;
    map[(int)'e'] = 9;
    map[(int)'k'] = 10;
    map[(int)'r'] = 11;
    map[(int)'n'] = 12;
    map[(int)'q'] = 13;
    map[(int)'c'] = 14;
    map[(int)'m'] = 15;
    map[(int)'w'] = 16;
    map[(int)'f'] = 17;
    map[(int)'y'] = 18;
    map[(int)'h'] = 19;
    map[(int)'.'] = 20;

    getInput(dnaSequence, proteinSequence, &dnaLen, &proteinLen);
    rdna = dnaLen - 4;
    rpro = proteinLen - 1;
    if (!findExon1(dnaSequence, proteinSequence, aminoacids, map, &ldna, &lpro))
    {
        printf("NONE");
        return 0;
    }
    intrn1Start = ldna;
    if (!findExon3(dnaSequence, proteinSequence, aminoacids, map, &rdna, &rpro))
    {
        printf("NONE");
        return 0;
    }
    intrn2End = rdna + 2;
    if(!(intrn1End = findExon2(dnaSequence, proteinSequence, aminoacids, map, ldna, lpro, rdna, rpro)-1)){
        printf("NONE");
        return 0;
    };
    intrn2Start = intrn1End + 1 + (rpro-lpro+1)*3;
    printf("\n%d %d %d %d", intrn1Start, intrn1End, intrn2Start, intrn2End);
    fflush(NULL);

    return 0;
}

int CheckCode(int key, char (*aminoacids)[6][3], char *readingFrame)
{
    int index = 0;
    while (aminoacids[key][index][0] != '\0')
    {
        if (aminoacids[key][index][0] == readingFrame[0] &&
            aminoacids[key][index][1] == readingFrame[1] &&
            aminoacids[key][index][2] == readingFrame[2])
        {
            return 1;
        }
        index++;
    }
    return 0;
}

int findExon1(char *dnaSequence, char *proteinSequence, char (*aminoacids)[6][3], int *map, int *ldna, int *lpro)
{
    char readingFrame[4] = "";
    while (dnaSequence[*ldna] != '\0')
    {
        if (dnaSequence[*ldna + 1] != '\0' && dnaSequence[*ldna + 2] != '\0')
        {
            readingFrame[0] = dnaSequence[*ldna];
            readingFrame[1] = dnaSequence[*ldna + 1];
            readingFrame[2] = dnaSequence[*ldna + 2];
            readingFrame[3] = '\0';
            if (CheckCode(map[(int)proteinSequence[*lpro]], aminoacids, readingFrame))
            {
                (*lpro)++;
            }
            else
            {
                return 1;
            }
        }
        (*ldna) += 3;
    }
    return 0;
}
int findExon2(char *dnaSequence, char *proteinSequence, char (*aminoacids)[6][3], int *map, int ldna, int lpro, int rdna, int rpro)
{
    char readingFrame[4] = "";
    int i;
    int k, l;
    for (; ldna < rdna; ldna++)
    {
        i = ldna;
        k = lpro;
        l = rpro;
        for (; k <= l; k++)
        {
            readingFrame[0] = dnaSequence[i];
            readingFrame[1] = dnaSequence[i + 1];
            readingFrame[2] = dnaSequence[i + 2];
            readingFrame[3] = '\0';
            if (!CheckCode(map[(int)proteinSequence[k]], aminoacids, readingFrame))
                break;
            i+=3;
        }
        if(k > l){
            return ldna;
        }
    }
    return 1;
}

int findExon3(char *dnaSequence, char *proteinSequence, char (*aminoacids)[6][3], int *map, int *rdna, int *rpro)
{
    char readingFrame[4] = "";
    readingFrame[0] = dnaSequence[*rdna];
    readingFrame[1] = dnaSequence[*rdna + 1];
    readingFrame[2] = dnaSequence[*rdna + 2];
    readingFrame[3] = '\0';
    if (!CheckCode(map[(int)proteinSequence[(*rpro)--]], aminoacids, readingFrame))
        return 0;
    for (; *rpro > 0; (*rpro)--)
    {
        (*rdna) -= 3;
        readingFrame[0] = dnaSequence[*rdna];
        readingFrame[1] = dnaSequence[*rdna + 1];
        readingFrame[2] = dnaSequence[*rdna + 2];
        readingFrame[3] = '\0';
        if (!CheckCode(map[(int)proteinSequence[*rpro]], aminoacids, readingFrame))
            return 1;
    }
    return 0;
}

void getInput(char *dnaSequence, char *proteinSequence, int *dnaLen, int *proteinLen)
{
    char ch;
    int count = 0;

    while ((ch = getchar()) != EOF)
    {
        if (ch == 'A' || ch == 'T' || ch == 'C' || ch == 'G')
        {
            dnaSequence[++(*dnaLen)] = ch;
        }
        else if (ch == 'g' || ch == 'a' || ch == 'v' || ch == 'i' || ch == 'l' || ch == 's' || ch == 't' || ch == 'p' || ch == 'd' || ch == 'e' || ch == 'k' || ch == 'r' || ch == 'n' || ch == 'q' || ch == 'c' || ch == 'm' || ch == 'w' || ch == 'f' || ch == 'y' || ch == 'h')
        {
            proteinSequence[++(*proteinLen)] = ch;
        }
        else if (ch == '.')
        {
            if (count == 0)
            {
                dnaSequence[++(*dnaLen)] = ch;
                count++;
            }
            else if (count == 1)
            {
                proteinSequence[++(*proteinLen)] = ch;
                count++;
            }
        }
        if (count == 2)
        {
            break;
        }
    }
    proteinSequence[++(*proteinLen)] = '\0';
    dnaSequence[++(*dnaLen)] = '\0';
}
