#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

static struct option long_options[] = {
    { "mode", required_argument, 0, 'm'},
    { "input", required_argument, 0, 'i'},
    { "gapopen", required_argument, 0, 'o'},
    { "gapext", required_argument, 0, 'e'}
};

enum directions { D, L, U};
enum modes { GLOBAL, AGLOBAL, LOCAL, ALOCAL};
char *outputFileNames[4] = { "global-naiveGap.aln", "global-affineGap.aln",
                            "local-naiveGap.aln", "local-affineGap.aln"};

int scoringMatrix[4][4] = {{2,-3,-3,-3},{-3,2,-3,-3},{-3,-3,2,-3},{-3,-3,-3,2}};
int openPen = -1;
int extPen = -1;

char **seqNames;
char **seqs;
int *seql;

int **scores;
enum directions **moves;

int alignmentLength;

int matchScore( int prev, char first, char second) {
    int fIndex = -1, sIndex = -1;

    switch( first) {
        case 'A':
            fIndex = 0;
            break;
        case 'C':
            fIndex = 1;
            break;
        case 'G':
            fIndex = 2;
            break;
        case 'T':
            fIndex = 3;
            break;
    }
    switch( second) {
        case 'A':
            sIndex = 0;
            break;
        case 'C':
            sIndex = 1;
            break;
        case 'G':
            sIndex = 2;
            break;
        case 'T':
            sIndex = 3;
            break;
    }

    return prev + scoringMatrix[fIndex][sIndex];
}

int naiveGapPenalty( int prev) {
    return prev + openPen;
}

int affineGapPenalty( int prev, int consecutive) {
    if ( consecutive == 1)
        return prev + extPen;
    return prev + openPen + extPen;
}

void globalAlignment() {
    scores[0][0] = 0;
    moves[0][0] = D;
    for ( int i = 1; i < seql[0]+1; i++) {
        scores[i][0] = naiveGapPenalty( scores[i-1][0]);
        moves[i][0] = U;
    }
    for ( int j = 1; j < seql[1]+1; j++) {
        scores[0][j] = naiveGapPenalty( scores[0][j-1]);
        moves[0][j] = L;
    }
    for ( int i = 1; i < seql[0]+1; i++) {
        for ( int j = 1; j < seql[1]+1; j++) {
            int maxScore, currentScore;
            enum directions move = D;

            maxScore = matchScore( scores[i-1][j-1], seqs[0][i-1], seqs[1][j-1]);

            currentScore = naiveGapPenalty( scores[i-1][j]);
            if ( currentScore > maxScore) {
                maxScore = currentScore;
                move = U;
            }

            currentScore = naiveGapPenalty( scores[i][j-1]);
            if ( currentScore > maxScore) {
                maxScore = currentScore;
                move = L;
            }

            scores[i][j] = maxScore;
            moves[i][j] = move;
        }
    }
}

void aGlobalAlignment() {
    scores[0][0] = 0;
    moves[0][0] = D;
    for ( int i = 1; i < seql[0]+1; i++) {
        if ( moves[i-1][0] == U)
            scores[i][0] = affineGapPenalty( scores[i-1][0], 1);
        else
            scores[i][0] = affineGapPenalty( scores[i-1][0], 0);
        moves[i][0] = U;
    }
    for ( int j = 1; j < seql[1]+1; j++) {
        if ( moves[0][j-1] == L)
            scores[0][j] = affineGapPenalty( scores[0][j-1], 1);
        else
            scores[0][j] = affineGapPenalty( scores[0][j-1], 0);
        moves[0][j] = L;
    }

    for ( int i = 1; i < seql[0]+1; i++) {
        for ( int j = 1; j < seql[1]+1; j++) {
            int maxScore, currentScore;
            enum directions move = D;
            maxScore = matchScore( scores[i-1][j-1], seqs[0][i-1], seqs[1][j-1]);

            int consecutive;
            if ( moves[i-1][j] == U)
                consecutive = 1;
            else
                consecutive = 0;
            currentScore = affineGapPenalty( scores[i-1][j], consecutive);
            if ( currentScore > maxScore) {
                maxScore = currentScore;
                move = U;
            }

            if ( moves[i][j-1] == L)
                consecutive = 1;
            else
                consecutive = 0;
            currentScore = affineGapPenalty( scores[i][j-1], consecutive);
            if ( currentScore > maxScore) {
                maxScore = currentScore;
                move = L;
            }
            scores[i][j] = maxScore;
            moves[i][j] = move;
        }
    }
}

void localAlignment() {
    scores[0][0] = 0;
    moves[0][0] = D;
    for ( int i = 1; i < seql[0]+1; i++) {
        scores[i][0] = naiveGapPenalty( scores[i-1][0]);
        moves[i][0] = U;
    }
    for ( int j = 1; j < seql[1]+1; j++) {
        scores[0][j] = naiveGapPenalty( scores[0][j-1]);
        moves[0][j] = L;
    }
    for ( int i = 1; i < seql[0]+1; i++) {
        for ( int j = 1; j < seql[1]+1; j++) {
            int maxScore, currentScore;
            enum directions move = D;

            maxScore = 0;

            currentScore = matchScore( scores[i-1][j-1], seqs[0][i-1], seqs[1][j-1]);
            if ( currentScore > maxScore)
                maxScore = currentScore;

            currentScore = naiveGapPenalty( scores[i-1][j]);
            if ( currentScore > maxScore) {
                maxScore = currentScore;
                move = U;
            }

            currentScore = naiveGapPenalty( scores[i][j-1]);
            if ( currentScore > maxScore) {
                maxScore = currentScore;
                move = L;
            }

            scores[i][j] = maxScore;
            moves[i][j] = move;
        }
    }
}

void aLocalAlignment() {
    scores[0][0] = 0;
    moves[0][0] = D;
    for ( int i = 1; i < seql[0]+1; i++) {
        if ( moves[i-1][0] == U)
            scores[i][0] = affineGapPenalty( scores[i-1][0], 1);
        else
            scores[i][0] = affineGapPenalty( scores[i-1][0], 0);
        moves[i][0] = U;
    }
    for ( int j = 1; j < seql[1]+1; j++) {
        if ( moves[0][j-1] == L)
            scores[0][j] = affineGapPenalty( scores[0][j-1], 1);
        else
            scores[0][j] = affineGapPenalty( scores[0][j-1], 0);
        moves[0][j] = L;
    }

    for ( int i = 1; i < seql[0]+1; i++) {
        for ( int j = 1; j < seql[1]+1; j++) {
            int maxScore, currentScore;
            enum directions move = D;

            maxScore = 0;

            currentScore = matchScore( scores[i-1][j-1], seqs[0][i-1], seqs[1][j-1]);
            if ( currentScore > maxScore)
                maxScore = currentScore;

            int consecutive;
            if ( moves[i-1][j] == U)
                consecutive = 1;
            else
                consecutive = 0;
            currentScore = affineGapPenalty( scores[i-1][j], consecutive);
            if ( currentScore > maxScore) {
                maxScore = currentScore;
                move = U;
            }

            if ( moves[i][j-1] == L)
                consecutive = 1;
            else
                consecutive = 0;
            currentScore = affineGapPenalty( scores[i][j-1], consecutive);
            if ( currentScore > maxScore) {
                maxScore = currentScore;
                move = L;
            }
            scores[i][j] = maxScore;
            moves[i][j] = move;
        }
    }
}

char *globalTraceback() {
    char *res, *temp;
    int index = 0;

    temp = (char *) malloc( seql[0] + seql[1] + 1);
    int i = seql[0], j = seql[1];
    while ( i >= 0 && j >= 0) {
        if ( moves[i][j] == D) {
            temp[index] = seqs[0][i-1];
            i--;
            j--;
        }
        else if ( moves[i][j] == U) {
            temp[index] = seqs[0][i-1];
            i--;
        }
        else if ( moves[i][j] == L) {
            temp[index] = seqs[1][j-1];
            j--;
        }
        index++;
    }
    index--;
    res = (char *) malloc( index);
    for ( int i = 0; i < index; i++) {
        res[index-1-i] = temp[i];
    }
    free(temp);

    alignmentLength = index;
    return res;
}

char *localTraceback( int x, int y) {
    char *res, *temp;
    int index = 0;

    temp = (char *) malloc( seql[0] + seql[1] + 1);
    int i = x, j = y;
    while ( scores[i][j] > 0 && i >= 0 && j >= 0) {
        if ( moves[i][j] == D) {
            temp[index] = seqs[0][i-1];
            i--;
            j--;
        }
        else if ( moves[i][j] == U) {
            temp[index] = seqs[0][i-1];
            i--;
        }
        else if ( moves[i][j] == L) {
            temp[index] = seqs[1][j-1];
            j--;
        }
        index++;
    }
    index--;
    res = (char *) malloc( index);
    for ( int i = 0; i < index; i++) {
        res[index-1-i] = temp[i];
    }
    free(temp);

    alignmentLength = index;
    return res;
}

char *threadSequenceGlobal( int n) {
    char *res, *temp;
    int index = 0;

    temp = (char *) malloc( seql[0] + seql[1] + 1);
    int i = seql[0], j = seql[1];
    while ( i >= 0 && j >= 0) {
        if ( moves[i][j] == D) {
            if ( n == 0)
                temp[index] = seqs[0][i-1];
            else
                temp[index] = seqs[1][j-1];
            i--;
            j--;
        }
        else if ( moves[i][j] == U) {
            if ( n == 0)
                temp[index] = seqs[0][i-1];
            else
                temp[index] = '-';
            i--;
        }
        else if ( moves[i][j] == L) {
            if ( n == 0)
                temp[index] = '-';
            else
                temp[index] = seqs[1][j-1];
            j--;
        }
        index++;
    }
    index--;
    res = (char *) malloc( index);
    for ( int i = 0; i < index; i++) {
        res[index-1-i] = temp[i];
    }
    free(temp);

    alignmentLength = index;
    return res;
}

char *threadSequenceLocal( int n, int x, int y) {
    char *res, *temp;
    int index = 0;

    temp = (char *) malloc( seql[0] + seql[1] + 1);
    int i = x, j = y;
    while ( scores[i][j] > 0 && i >= 0 && j >= 0) {
        if ( moves[i][j] == D) {
            if ( n == 0)
                temp[index] = seqs[0][i-1];
            else
                temp[index] = seqs[1][j-1];
            i--;
            j--;
        }
        else if ( moves[i][j] == U) {
            if ( n == 0)
                temp[index] = seqs[0][i-1];
            else
                temp[index] = '-';
            i--;
        }
        else if ( moves[i][j] == L) {
            if ( n == 0)
                temp[index] = '-';
            else
                temp[index] = seqs[1][j-1];
            j--;
        }
        index++;
    }
    index--;
    res = (char *) malloc( index);
    for ( int i = 0; i < index; i++) {
        res[index-1-i] = temp[i];
    }
    free(temp);

    alignmentLength = index;
    return res;
}

long ms( struct timeval t) {
    return t.tv_sec * 1000000 + t.tv_usec;
}

int main( int argc, char** argv)
{
    char *inputFile;

    enum modes mode;
    int index;
    int opt;
    int mFlag = 0, iFlag = 0;

    while ((opt = getopt_long( argc, argv, "m:i:o:e", long_options, &index)) != -1) {
        switch(opt) {
            case 'm':
                if ( strcmp( "global", optarg) == 0)
                    mode = GLOBAL;
                else if  ( strcmp( "local", optarg) == 0)
                    mode = LOCAL;
                else if ( strcmp( "aglobal", optarg) == 0)
                    mode = AGLOBAL;
                else if  ( strcmp( "alocal", optarg) == 0)
                    mode = ALOCAL;
                mFlag = 1;
                break;
            case 'i':
                inputFile = (char *) malloc( 25);
                strcpy( inputFile, optarg);
                iFlag = 1;
                break;
            case 'o':
                openPen = atoi( optarg);
                break;
            case 'e':
                extPen = atoi( optarg);
                break;
        }
    }

    if ( iFlag == 0) {
        printf( "Input file missing.\n");
        return 1;
    }

    if ( mFlag == 0) {
        printf( "Mode selection missing, global mode selected.\n");
        mode = GLOBAL;
    }

    // allocate memory
    seqNames = (char **) malloc( 2 * sizeof(char *));
    seqs = (char **) malloc( 2 * sizeof(char *));
    seql = (int *) malloc( 2 * sizeof(int));

    // read input file
    int seqIndex = -1;
    char * line;
    size_t len = 0;
    FILE *fp = fopen( inputFile, "r");
    if ( fp == NULL) {
        printf( "Can't open \"%s|"".\n", inputFile);
        return 1;
    }
    while ((getline(&line, &len, fp)) != -1) {
        if ( line[0] == 62) {
            seqIndex++;

            int nameLen = 0;
            for ( int i = 0; i < len; i++)
                if ( line[i] != 10)
                    nameLen++;
                else
                    break;
            seqNames[seqIndex] = (char *) malloc( nameLen);
            for ( int i = 0; i < nameLen-1; i++)
            {
                seqNames[seqIndex][i] = line[i+1];
            }
        }
        else {
            seql[seqIndex] = 0;
            for ( int i = 0; i < len; i++)
                if ( line[i] == 65
                    || line[i] == 67
                    || line[i] == 71
                    || line[i] == 84)
                    seql[seqIndex]++;
                else
                    break;
            seqs[seqIndex] = malloc( seql[seqIndex] + 1);
            for ( int i = 0; i < seql[seqIndex]; i++)
                seqs[seqIndex][i] = line[i];
        }
    }
    fclose(fp);
    if (line)
        free(line);
    free(inputFile);

    // create allignment matrix
    scores = (int **) malloc( (1 + seql[0]) * sizeof(int *));
    for ( int i = 0; i < seql[0]+1; i++) {
        scores[i] = (int *) malloc( (1 + seql[1]) * sizeof(int));
    }
    moves = (enum directions **) malloc( (1 + seql[0]) * sizeof(enum directions *));
    for ( int i = 0; i < seql[0]+1; i++) {
        moves[i] = (enum directions *) malloc( (1 + seql[1]) * sizeof(enum directions));
    }

    // calculate scores according to mode
    int oFileIndex = 0;
    switch (mode) {
        case GLOBAL:
            globalAlignment();
            oFileIndex = 0;
            break;
        case AGLOBAL:
            aGlobalAlignment();
            oFileIndex = 1;
            break;
        case LOCAL:
            localAlignment();
            oFileIndex = 2;
            break;
        case ALOCAL:
            aLocalAlignment();
            oFileIndex = 3;
            break;
    }

    char **alignedSeqs = (char **) malloc( 2 * sizeof(char *));

    if ( mode == GLOBAL || mode == AGLOBAL) {
        alignedSeqs[0] = threadSequenceGlobal(0);
        alignedSeqs[1] = threadSequenceGlobal(1);
    }
    else {
        int x = seql[0], y = seql[1];
        int maxScore = scores[seql[0]][seql[1]];
        for ( int i = seql[0]; i >= 0; i--) {
            for ( int j = seql[1]; j  >= 0; j--) {
                if ( scores[i][j] > maxScore) {
                    maxScore = scores[i][j];
                    x = i;
                    y = j;
                }
            }
        }
        alignedSeqs[0] = threadSequenceLocal( 0, x, y);
        alignedSeqs[1] = threadSequenceLocal( 1, x, y);
    }


    // open oputput file
    fp = fopen( outputFileNames[oFileIndex], "w");
    if ( fp == NULL) {
        printf( "Unable to open output file.\n");
        return 1;
    }

    // write results
    fprintf( fp, "Score = %d\n", scores[seql[0]][seql[1]]);

    int j = 0, n = 0;
    for ( int i = 0; i <= alignmentLength; i++) {
        if ( i == alignmentLength && n > 0) {
            break;
        }
        if ( j >= 60 || i == alignmentLength) {
            fprintf( fp, "\n");
            if ( n == 0) {
                i = i - j;
                n++;
            }
            else
                n = 0;
            j = 0;
        }
        if ( j == 0)
            fprintf( fp, "%-20s ", seqNames[n]);
        j++;
        fprintf( fp, "%c", alignedSeqs[n][i]);
    }
    fclose( fp);

    // memory deallocation
    for ( int i = 0; i < 2; i++) {
        free( seqNames[i]);
        free( seqs[i]);
        free( alignedSeqs[i]);
    }
    free( seqNames);
    free( seqs);
    free( alignedSeqs);
    for ( int i = 0; i < seql[0]+1; i++) {
        free(scores[i]);
        free(moves[i]);
    }
    free( scores);
    free( moves);
    free( seql);


/*
    struct timeval t1, t2;
    int bFComp, bFPos;
    long bFTime;

    gettimeofday( &t1, NULL);


    gettimeofday( &t2, NULL);
    bFTime = (ms(t2) - ms(t1));

    printf( "Runtime was %ldms.\n", bFTime);
*/

    return 0;
}
