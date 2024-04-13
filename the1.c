#include <stdio.h>
#define MAX_DNA_LENGTH 3005
#define MAX_PROTEIN_LENGTH 3000

const char aminoacids[21][6][3] = {
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
void getInput(char *, char *, int *, int *);
int CheckCode(int key, char *readingFrame);
int findExon2(char *, char *, int *, int, int, int, int);
int findExon3(char *, char *, int *, int *, int *);

int main(void)
{
    int map[200];
    int intrn1Start, intrn2Start, intrn1End, intrn2End;
    char dnaSequence[MAX_DNA_LENGTH] = "";
    char proteinSequence[MAX_PROTEIN_LENGTH] = "";
    int proteinLen = -1;
    int dnaLen = -1;
    int ldna, rdna, lpro, rpro;
    char readingFrame[4];
    int check = 1;
    int k, l;
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
    rdna = dnaLen - 3;
    rpro = proteinLen - 1;

    /*Checking stop codon*/
    readingFrame[0] = dnaSequence[rdna];
    readingFrame[1] = dnaSequence[rdna + 1];
    readingFrame[2] = dnaSequence[rdna + 2];
    readingFrame[3] = '\0';
    if (!CheckCode(map['.'], readingFrame))
    {
        printf("NONE");
        return 0;
    }

    if (findExon3(dnaSequence, proteinSequence, map, &rdna, &rpro) == -1)
    {
        if(rpro < 0)
        {
            rpro += 2;
            rdna += 5;
        }
    }
    intrn2End = rdna;

    /*Check First Codon*/
    readingFrame[0] = dnaSequence[ldna];
    readingFrame[1] = dnaSequence[ldna + 1];
    readingFrame[2] = dnaSequence[ldna + 2];
    readingFrame[3] = '\0';
    if (CheckCode(map[(int)proteinSequence[lpro]], readingFrame))
    {
        ldna += 3;
        lpro++;
    }
    else
    {
        printf("NONE");
        return 0;
    }

    while (1)
    {
        readingFrame[0] = dnaSequence[ldna];
        readingFrame[1] = dnaSequence[ldna + 1];
        readingFrame[2] = dnaSequence[ldna + 2];
        readingFrame[3] = '\0';
        if (CheckCode(map[(int)proteinSequence[lpro]], readingFrame) || check)
        {
            if (!CheckCode(map[(int)proteinSequence[lpro]], readingFrame))
            {
                check = 0;
            }
            if (findExon2(dnaSequence, proteinSequence, map, ldna + 1, lpro, rdna, rpro) != -1)
            {
                intrn1Start = ldna;
                intrn1End = map[1];
                intrn2Start = map[0];
                printf("%d %d %d %d", intrn1Start, intrn1End, intrn2Start, intrn2End);
                fflush(NULL);
                if (intrn1Start > intrn1End || intrn2Start > intrn2End || intrn1Start > intrn2Start || intrn1End > intrn2End)
                {
                    printf("NONE");
                    fflush(NULL);
                    return 0;
                }
                return 0;
            }
            else
            {
                k = rpro;
                l = rdna;
                while (k < proteinLen-1)
                {
                    k++;
                    l += 3;
                    if (findExon2(dnaSequence, proteinSequence, map, ldna + 1, lpro, l, k) != -1)
                    {
                        intrn1Start = ldna;
                        intrn1End = map[1];
                        intrn2Start = map[0];
                        intrn2End = l;
                        if (intrn1Start > intrn1End || intrn2Start > intrn2End || intrn1Start > intrn2Start || intrn1End > intrn2End)
                        {
                            printf("NONE");
                            fflush(NULL);
                            return 0;
                        }
                        printf("%d %d %d %d", intrn1Start, intrn1End, intrn2Start, intrn2End);
                        fflush(NULL);
                        return 0;
                    }
                }
            }
            lpro++;
            ldna += 3;
        }
        else
        {
            printf("NONE");
            fflush(NULL);
            return 0;
        }
    }
}

int findExon2(char *dnaSequence, char *proteinSequence, int *map, int ldna, int lpro, int rdna, int rpro)
{
    char readingFrame[4] = "";
    int i;
    int k, l;
    for (; ldna + 3 * (rpro - lpro) < rdna; ldna++)
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
            if (!CheckCode(map[(int)proteinSequence[k]], readingFrame))
                break;
            i += 3;
        }
        if (k > l && lpro <= rpro)
        {
            map[0] = i;
            map[1] = ldna - 1;
            return ldna - 1;
        }
    }
    return -1;
}

int CheckCode(int key, char *readingFrame)
{
    int index = 0;
    while (aminoacids[key][index][0] != '\0' && index < 6)
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

int findExon3(char *dnaSequence, char *proteinSequence, int *map, int *rdna, int *rpro)
{
    char readingFrame[4] = "";
    for (; *rpro >= 0; (*rpro)--)
    {
        (*rdna) -= 3;
        readingFrame[0] = dnaSequence[*rdna];
        readingFrame[1] = dnaSequence[*rdna + 1];
        readingFrame[2] = dnaSequence[*rdna + 2];
        readingFrame[3] = '\0';
        if (!CheckCode(map[(int)proteinSequence[*rpro]], readingFrame))
        {
            (*rdna) += 2;
            return 1;
        }
    }
    return -1;
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
                count++;
            }
            else if (count == 1)
            {
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