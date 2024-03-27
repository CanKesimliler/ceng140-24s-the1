#include <stdio.h>
#define MAX_DNA_LENGTH 3005
#define MAX_PROTEIN_LENGTH 3000 

int CheckCode(int key, char*** aminoacids, char* readingFrame);

int main(){

    short int map[100] = {};
    map['g'] = 0;
    map['a'] = 1;
    map['v'] = 2;
    map['i'] = 3;
    map['l'] = 4;
    map['s'] = 5;
    map['t'] = 6;
    map['p'] = 7;
    map['d'] = 8;
    map['e'] = 9;
    map['k'] = 10;
    map['r'] = 11;
    map['n'] = 12;
    map['q'] = 13;
    map['c'] = 14;
    map['m'] = 15;
    map['w'] = 16;
    map['f'] = 17;
    map['y'] = 18;
    map['h'] = 19;
    map['.'] = 20;
    char aminoacids[21][6][3] = {
        {"GGT", "GGC", "GGA", "GGG", "", ""},
        {"GCT", "GCC", "GCA", "GCG", "", ""},
        {"GTT", "GTC", "GTA", "GTG", "", ""},
        {"ATT", "ATC", "ATA", "", "", ""},
        {"TTA", "TTG" ,"CTT", "CTC", "CTA", "CTG"},
        {"AGT", "AGC", "TCT", "TCC", "TCA", "TCG"},
        {"ACT", "ACC", "ACA", "ACG", "", ""},
        {"CCT", "CCC", "CCA", "CCG", "", ""},
        {"GAT", "GAG", "", "", "", ""},
        {"GAA", "GAG", "", "", "", ""},
        {"AAA", "AAG", "", "", "", ""},
        {"AGA", "AGG", "CGT", "CGC", "CGA", "CGG"},
        {"AAT", "AAC", "", "", "", ""},
        {"CAA", "CAG", "", "", "", ""},
        {"TGT", "TGC", "", "", "", ""},
        {"ATG", "", "", "", "", ""},
        {"TGG", "", "", "", "", ""},
        {"TTT", "TTC", "", "", "", ""},
        {"TAT", "TAC", "", "", "", ""},
        {"CAT", "CAC", "", "", "", ""},
        {"TAA", "TAG", "TGA", "", "", ""}
    };
    int i = 0;
    int k = 0;
    int intrn1Start, intrn2Start, intrn1End, intrn2End;
    intrn1Start = intrn1End = intrn2Start = intrn2End = 0;
    char dnaSequence[MAX_DNA_LENGTH];
    char proteinSequence[MAX_PROTEIN_LENGTH];
    char readingFrame[4];

    scanf("%s %s",&dnaSequence, &proteinSequence);

    while(dnaSequence[i] != '\0'){
        if(dnaSequence[i+1] != '\0' && dnaSequence[i+2] != '\0'){
            readingFrame[0] = dnaSequence[i];
            readingFrame[1] = dnaSequence[i+1];
            readingFrame[2] = dnaSequence[i+2];
            readingFrame[3] = '\0';
            if(CheckCode(map[proteinSequence[k]], aminoacids, readingFrame)){
                k++;
            }
            else if(intrn1Start == 0){
                
            }
        }
    }




    return 0;
}

int CheckCode(int key, char*** aminoacids, char* readingFrame){
    int index = 0;
    int i;
    while(aminoacids[key][index][0] != '\0'){
        if(aminoacids[key][index][0] == readingFrame[0] && 
        aminoacids[key][index][1] == readingFrame[1] && 
        aminoacids[key][index][2] == readingFrame[2]){
            return 1;
        }
        index++;
    }
    return 0;
}