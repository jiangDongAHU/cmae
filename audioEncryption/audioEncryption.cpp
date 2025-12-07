#include "audioEncryption.hpp"

int main(){
    
    //configuration structure for initializing capture and palyback devices
    static const pa_sample_spec deviceConfigurationStructure = {
        .format   = PA_SAMPLE_S16LE, // 16-bit signed, little-endian
        .rate     = SAMPLING_RATE,   // sampling rate in Hz
        .channels = CHANNELS         // number of channels
    };

    pa_simple * inputStream  = NULL;
    pa_simple * outputStream = NULL;
    int         err;

    //create the input stream
    if(!(inputStream = pa_simple_new(NULL, "audioEncryption", PA_STREAM_RECORD, NULL, "capture", &deviceConfigurationStructure, NULL, NULL, NULL))){
        fprintf(stderr, "unable to create the input stream: %s\n", pa_strerror(err));
        exit(1);
    }

    //create the output stream
    if(!(outputStream = pa_simple_new(NULL, "audioEncryption", PA_STREAM_PLAYBACK, NULL, "playback", &deviceConfigurationStructure, NULL, NULL, NULL))){
        fprintf(stderr, "unable to create the output stream: %s\n", pa_strerror(err));
        exit(1);
    }

    //initialize SDL
    if(SDL_Init(SDL_INIT_VIDEO) < 0){
        fprintf(stderr, "unable to initialize SDL: %s\n", SDL_GetError());
        exit(1);
    }

    SDL_Window *window = SDL_CreateWindow(
        "Original, encrypted, and decrypted audio samples",
        SDL_WINDOWPOS_UNDEFINED,
        SDL_WINDOWPOS_UNDEFINED,
        WINDOW_WIDTH,
        WINDOW_HEIGHT,
        SDL_WINDOW_SHOWN
    );

    if (!window) {
        fprintf(stderr, "unable to create window: %s\n", SDL_GetError());
        return 1;
    }

    SDL_Renderer * renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (!renderer) {
        fprintf(stderr, "unable to create renderer: %s\n", SDL_GetError());
        return 1;
    }

    //randomly select initial conditions and control parameters as the user inputs
    srand(time(NULL));
    double initialCondition1 = (double)rand() / RAND_MAX;
    double controlParameter1 = (double)rand() / RAND_MAX;
    double initialCondition2 = (double)rand() / RAND_MAX;
    double controlParameter2 = (double)rand() / RAND_MAX;
    if(controlParameter1 > 0.5)
        controlParameter1 = 1 - controlParameter1;
    if(controlParameter2 > 0.5)
        controlParameter2 = 1 - controlParameter2;

    //allocate memory for processing audio samples
    uint8_t * sampleBuffer = (uint8_t *)malloc(SAMPLE_BUFFER_SIZE * sizeof(uint8_t));
    uint8_t SHA256HashResultsArray[SHA256_DIGEST_LENGTH];
    uint8_t originalByteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS];
    uint8_t tempByteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS];
    uint8_t encryptedByteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS];
    uint8_t decryptedByteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS];

    int iterationsForGeneratingShiftDistanceSequence = (int)((BYTE_MATRIX_ROWS + BYTE_MATRIX_COLUMNS) * 
                                                             (ROUNDS_OF_CONFUSION_OPERATIONS + ROUNDS_OF_DIFFUSION_CONFUSION_OPERATIONS) /
                                                              NUMBER_OF_UINT16T_DATA_EXTRACTED_FROM_EACH_ITERATION_RESULT) + 1;
    double * iterationResultsArray1ForGeneratingShiftDistanceSequence = (double *)malloc(iterationsForGeneratingShiftDistanceSequence * sizeof(double));
    double * iterationResultsArray2ForGeneratingShiftDistanceSequence = (double *)malloc(iterationsForGeneratingShiftDistanceSequence * sizeof(double));
    uint16_t * shiftDistanceSequence = (uint16_t *)malloc(iterationsForGeneratingShiftDistanceSequence * NUMBER_OF_UINT16T_DATA_EXTRACTED_FROM_EACH_ITERATION_RESULT * sizeof(uint16_t));

    int iterationsForGeneratingByteSequence = (int)((BYTE_MATRIX_ROWS * BYTE_MATRIX_COLUMNS) * ROUNDS_OF_DIFFUSION_CONFUSION_OPERATIONS / 
                                                     NUMBER_OF_UINT8T_DATA_EXTRACTED_FROM_EACH_ITERATION_RESULT) + 1;
    double * iterationResultArray1ForGeneratingByteSequence = (double *)malloc(iterationsForGeneratingByteSequence * sizeof(double));
    double * iterationResultArray2ForGeneratingByteSequence = (double *)malloc(iterationsForGeneratingByteSequence * sizeof(double));
    uint8_t * byteSequence = (uint8_t *)malloc(iterationsForGeneratingByteSequence * NUMBER_OF_UINT8T_DATA_EXTRACTED_FROM_EACH_ITERATION_RESULT * sizeof(uint8_t));

    Mat originalFrame  = Mat::zeros(BYTE_MATRIX_ROWS, BYTE_MATRIX_COLUMNS, CV_8UC1);
    Mat encryptedFrame = Mat::zeros(BYTE_MATRIX_ROWS, BYTE_MATRIX_COLUMNS, CV_8UC1);
    Mat decryptedFrame = Mat::zeros(BYTE_MATRIX_ROWS, BYTE_MATRIX_COLUMNS, CV_8UC1);
    Mat canvas         = Mat::zeros(WINDOW_HEIGHT,    WINDOW_WIDTH,        CV_8UC1);

    int blockIndex = 0;
    while(blockIndex < 10000){
        //capture samples
        if(pa_simple_read(inputStream, sampleBuffer, SAMPLE_BUFFER_SIZE, &err) < 0){
            fprintf(stderr, "unable to read samples: %s\n", pa_strerror(err));
            exit(1);
        }

        double startTime = getCPUSecond();

        //calculate the SHA-256 hash of the captured samples
        calculateSHA256Hash(sampleBuffer, SHA256HashResultsArray);
        initialCondition1 = reconstructParameter(initialCondition1, &SHA256HashResultsArray[0]);
        controlParameter1 = reconstructParameter(controlParameter1, &SHA256HashResultsArray[8]);
        initialCondition2 = reconstructParameter(initialCondition2, &SHA256HashResultsArray[16]);
        controlParameter2 = reconstructParameter(controlParameter2, &SHA256HashResultsArray[24]);
        if(controlParameter1 > 0.5)
            controlParameter1 = 1 - controlParameter1;
        if(controlParameter2 > 0.5)
            controlParameter2 = 1 - controlParameter2;

        //pre-iterate PLCMs
        for(int i = 0; i < PRE_ITERATIONS; i++){
            initialCondition1 = PLCM(initialCondition1, controlParameter1);
            initialCondition2 = PLCM(initialCondition2, controlParameter2);
        }

        //generate iteration results for generating shif distance sequence
        for(int i = 0; i < iterationsForGeneratingShiftDistanceSequence; i++){
            initialCondition1 = PLCM(initialCondition1, controlParameter1);
            iterationResultsArray1ForGeneratingShiftDistanceSequence[i] = initialCondition1;

            initialCondition2 = PLCM(initialCondition2, controlParameter2);
            iterationResultsArray2ForGeneratingShiftDistanceSequence[i] = initialCondition2;
        }

        //construct shift distance sequence
        convertIterationResultsToShiftDistanceSequence(iterationResultsArray1ForGeneratingShiftDistanceSequence,
                                                       iterationResultsArray2ForGeneratingShiftDistanceSequence,
                                                       shiftDistanceSequence, 
                                                       iterationsForGeneratingShiftDistanceSequence);
        
        //generate iteration results for generating byte sequence
        for(int i = 0; i < iterationsForGeneratingByteSequence; i++){
            initialCondition1 = PLCM(initialCondition1, controlParameter1);
            iterationResultArray1ForGeneratingByteSequence[i] = initialCondition1;

            initialCondition2 = PLCM(initialCondition2, controlParameter2);
            iterationResultArray2ForGeneratingByteSequence[i] = initialCondition2;
        }

        //construct byte sequence
        convertIterationResultsToByteSequence(iterationResultArray1ForGeneratingByteSequence,
                                              iterationResultArray2ForGeneratingByteSequence,
                                              byteSequence,
                                              iterationsForGeneratingByteSequence);
 
        //convert the linear samples into byte matrix
        convertSampleBufferToByteMatrix(sampleBuffer, originalByteMatrix);
        cloneByteMatrix(originalByteMatrix, tempByteMatrix);

        //perform confusion operations on the byte matrix
        int shiftDistanceSequenceIndex = 0;
        for(int i = 0; i < ROUNDS_OF_CONFUSION_OPERATIONS; i++)
            shiftDistanceSequenceIndex = confusion(encryptedByteMatrix, tempByteMatrix, shiftDistanceSequence, shiftDistanceSequenceIndex);

        //perform diffusion and confusion operations on the byte matrix
        int byteSequenceIndex = 0;
        uint8_t diffusionSeedArray[ROUNDS_OF_DIFFUSION_CONFUSION_OPERATIONS];

        for(int i = 0; i < ROUNDS_OF_CONFUSION_OPERATIONS; i++){
            //perform diffusion operations
            diffusionSeedArray[i] = generateDiffusionSeed(&initialCondition1, controlParameter1, &initialCondition2, controlParameter2);
            byteSequenceIndex     = diffusion(encryptedByteMatrix, tempByteMatrix, byteSequence, diffusionSeedArray[i], byteSequenceIndex);

            //perform confusion operations
            shiftDistanceSequenceIndex = confusion(encryptedByteMatrix, tempByteMatrix, shiftDistanceSequence, shiftDistanceSequenceIndex);
        }

        double encryptionTime = getCPUSecond();

        //perform inverse confusion and diffusion operations on the byte matrix
        for(int i = ROUNDS_OF_DIFFUSION_CONFUSION_OPERATIONS - 1; i >= 0; i--){
            //perform inverse confusion operations
            shiftDistanceSequenceIndex = inverseConfusion(decryptedByteMatrix, tempByteMatrix, shiftDistanceSequence, shiftDistanceSequenceIndex);

            //perform inverse diffusion operations
            byteSequenceIndex = inverseDiffusion(decryptedByteMatrix, tempByteMatrix, byteSequence, diffusionSeedArray[i], byteSequenceIndex);
        }

        //perform inverse confusion operations on the byte matrix
        for(int i = 0; i < ROUNDS_OF_CONFUSION_OPERATIONS; i++)
            shiftDistanceSequenceIndex = inverseConfusion(decryptedByteMatrix, tempByteMatrix, shiftDistanceSequence, shiftDistanceSequenceIndex);

        double endTime = getCPUSecond();

        //draw the waveform of the audio samples
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255); 
        SDL_RenderClear(renderer);

        drawWaveform(renderer, sampleBuffer, SAMPLE_BUFFER_SIZE, 0);

        //convert the byte matrix into linear samples
        convertByteMatrixToSampleBuffer(encryptedByteMatrix, sampleBuffer);
        drawWaveform(renderer, sampleBuffer, SAMPLE_BUFFER_SIZE, 1);

        convertByteMatrixToSampleBuffer(decryptedByteMatrix, sampleBuffer);
        drawWaveform(renderer, sampleBuffer, SAMPLE_BUFFER_SIZE, 2);
        SDL_RenderPresent(renderer);

        //draw the byte matrix
        for(int i = 0; i < BYTE_MATRIX_ROWS; i++)
            for(int j = 0; j < BYTE_MATRIX_COLUMNS; j++){
                originalFrame.at<uchar>(i,j)  = originalByteMatrix[i][j];
                encryptedFrame.at<uchar>(i,j) = encryptedByteMatrix[i][j];
                decryptedFrame.at<uchar>(i,j) = decryptedByteMatrix[i][j];
            }

        drawByteMatrix(originalFrame, canvas, 0);
        drawByteMatrix(encryptedFrame, canvas, 1);
        drawByteMatrix(decryptedFrame, canvas, 2);

        imshow("original, encrypted, and decrypted byte matrix", canvas);
        waitKey(1);

        //playback the captured samples
        if(pa_simple_write(outputStream, sampleBuffer, SAMPLE_BUFFER_SIZE, &err) < 0){
            fprintf(stderr, "unable to playback samples: %s\n", pa_strerror(err));
            exit(1);
        }

              
        if(blockIndex % 50 == 0){
            system("clear");
            printf("Device and Protocol Configurations:\n");
            printf("Channels: %d  || Sample Width: %d bytes   || Sampling Rate: %d Hz\n", CHANNELS, SAMPLE_WIDTH, SAMPLING_RATE);
            printf("Frames: %d || Rows of Byte Matrix: %d || Columns of Byte Matrix: %d\n", FRAMES_IN_EACH_BLOCK, BYTE_MATRIX_ROWS, BYTE_MATRIX_COLUMNS);
            printf("Rounds of Confusion Operations: %d       || Rounds of Diffusion and Confusion Operations: %d\n\n", ROUNDS_OF_CONFUSION_OPERATIONS, ROUNDS_OF_DIFFUSION_CONFUSION_OPERATIONS);
            
            printf("SHA-256 Hash of the Original Samples:\n");
            dispalySHA256Hash(SHA256HashResultsArray);

            printf("\nReconstructed Initial Conditions and Control Parameters:\n");
            printf("Initial Condition 1: %.6f           || Control Parameter 1: %.6f\n", initialCondition1, controlParameter1);
            printf("Initial Condition 2: %.6f           || Control Parameter 2: %.6f\n\n", initialCondition2, controlParameter2);

            printf("Encryption and Decryption Speed:\n");
            printf("Encryption Time: %.2f ms                || Encryption and Decryption Time: %.2f ms  \n", (encryptionTime - startTime) * 1000, (endTime - startTime) * 1000);
        }

        blockIndex++;
    }

    free(sampleBuffer);
    free(iterationResultsArray1ForGeneratingShiftDistanceSequence);
    free(iterationResultsArray2ForGeneratingShiftDistanceSequence);
    free(shiftDistanceSequence);
    free(iterationResultArray1ForGeneratingByteSequence);
    free(iterationResultArray2ForGeneratingByteSequence);
    free(byteSequence);
    
    return 0;
}

//get cpu time
double getCPUSecond(void){
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}


//convert the linear samples into byte matrix
void convertSampleBufferToByteMatrix(uint8_t * sampleBuffer, uint8_t byteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS]){
    int sampleBufferIndex = 0;

    for(int i = 0; i < BYTE_MATRIX_ROWS; i++)
        for(int j = 0; j < BYTE_MATRIX_COLUMNS; j++)
            byteMatrix[i][j] = sampleBuffer[sampleBufferIndex ++];
}

//convert the byte matrix into the linear samples
void convertByteMatrixToSampleBuffer(uint8_t byteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS], uint8_t * sampleBuffer){
    int sampleBufferIndex = 0;

    for(int i = 0; i < BYTE_MATRIX_ROWS; i++)
        for(int j = 0; j < BYTE_MATRIX_COLUMNS ;j++)
        sampleBuffer[sampleBufferIndex ++] = byteMatrix[i][j];
}

//calculate the SHA-256 hash
void calculateSHA256Hash(uint8_t * sampleBuffer, uint8_t * SHA256HashResultsArray){
    //SHA-256 structure
    SHA256_CTX sha256;

    //initialize the sturcture
    SHA256_Init(& sha256);

    //update the structure with the samples
    SHA256_Update(& sha256, sampleBuffer, SAMPLE_BUFFER_SIZE);

    //finalize the hash and store the result
    SHA256_Final(SHA256HashResultsArray, & sha256);
}

//reconsturct initial condition or control parameter using SHA-256 hash of the captured samples
double reconstructParameter(double parameter, uint8_t * SHA256HashResultsArray){
    uint64_t uint64Parameter = 0;
    memcpy(&uint64Parameter, &parameter, sizeof(uint64_t));

    //fetch 8 bytes from the SHA256 hash
    uint64_t SHA256HashValue = 0;
    memcpy(&SHA256HashValue, SHA256HashResultsArray, 8);

    //perform XOR operations
    uint64_t data = uint64Parameter ^ SHA256HashValue;

    double result = (double)data / 0x1p64;

    return result;
}

//iterate the PLCM and return the iteration result
double PLCM(double initialCondition, double controlParameter){
    double iterationResult = 0;

    if(initialCondition >= 0 && initialCondition <= controlParameter)
        iterationResult = initialCondition / controlParameter;
    
    else if(initialCondition > controlParameter && initialCondition <= 0.5)
        iterationResult = (initialCondition - controlParameter) / (0.5 - controlParameter);
    
    else
        iterationResult = PLCM(controlParameter, 1 - initialCondition);

    return iterationResult;
}

//convert iteration results into shift distance sequence
void convertIterationResultsToShiftDistanceSequence(double * iterationResultArray1, double * iterationResultArray2, uint16_t * shiftDistanceSequence, int numberOfIterationResults){
    uint16_t uint16tDataArray1[numberOfIterationResults * NUMBER_OF_UINT16T_DATA_EXTRACTED_FROM_EACH_ITERATION_RESULT];
    uint16_t uint16tDataArray2[numberOfIterationResults * NUMBER_OF_UINT16T_DATA_EXTRACTED_FROM_EACH_ITERATION_RESULT];

    //convert iteration results into uint16_t data
    for(int i = 0; i < numberOfIterationResults; i++){
        uint16_t * p = &uint16tDataArray1[i * NUMBER_OF_UINT16T_DATA_EXTRACTED_FROM_EACH_ITERATION_RESULT];
        memcpy(p, &iterationResultArray1[i], NUMBER_OF_UINT16T_DATA_EXTRACTED_FROM_EACH_ITERATION_RESULT);
        
        p = &uint16tDataArray2[i * NUMBER_OF_UINT16T_DATA_EXTRACTED_FROM_EACH_ITERATION_RESULT];
        memcpy(p, &iterationResultArray2[i], NUMBER_OF_UINT16T_DATA_EXTRACTED_FROM_EACH_ITERATION_RESULT);
    }

    //generate shift distances
    for(int i = 0; i < numberOfIterationResults * NUMBER_OF_UINT16T_DATA_EXTRACTED_FROM_EACH_ITERATION_RESULT; i++){
        shiftDistanceSequence[i] = uint16tDataArray1[i] ^ uint16tDataArray2[i];
    }
}

//clone byte matrix
void cloneByteMatrix(uint8_t sourceByteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS], uint8_t targetBytematrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS]){
    for(int i = 0; i < BYTE_MATRIX_ROWS; i++)
        for(int j = 0; j < BYTE_MATRIX_COLUMNS; j++)
            targetBytematrix[i][j] = sourceByteMatrix[i][j];
}

//perform confusion operations on the byte matrix
int confusion(uint8_t byteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS], uint8_t tempByteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS], uint16_t * shiftDistanceSequence, int shiftDistanceSequenceIndex){   
    uint16_t byteMatrixColumns = BYTE_MATRIX_COLUMNS;
    uint16_t byteMatrixRows    = BYTE_MATRIX_ROWS;

    //perform confusion operations along horizontal direction
    for(int i = 0; i < BYTE_MATRIX_ROWS; i++){
        uint16_t shiftDistance = shiftDistanceSequence[shiftDistanceSequenceIndex ++] % byteMatrixColumns; 

        for(int j = 0; j < BYTE_MATRIX_COLUMNS; j++)
            byteMatrix[i][(j + shiftDistance) % byteMatrixColumns] = tempByteMatrix[i][j];        
    }

    cloneByteMatrix(byteMatrix, tempByteMatrix);

    //perform confusion operations along vertical direction
    for(int j = 0; j < BYTE_MATRIX_COLUMNS; j++){
        uint16_t shiftDistance = shiftDistanceSequence[shiftDistanceSequenceIndex ++] % byteMatrixRows; 

        for(int i = 0; i < BYTE_MATRIX_ROWS; i++)
            byteMatrix[(i + shiftDistance) % byteMatrixRows][j] = tempByteMatrix[i][j];
    }
    cloneByteMatrix(byteMatrix, tempByteMatrix);

    return shiftDistanceSequenceIndex;
}

//perform inverse confusion operations on the byte matrix
int inverseConfusion(uint8_t byteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS], uint8_t tempByteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS], uint16_t * shiftDistanceSequence, int shiftDistanceSequenceIndex){
    uint16_t byteMatrixColumns = BYTE_MATRIX_COLUMNS;
    uint16_t byteMatrixRows    = BYTE_MATRIX_ROWS;

    //perform inverse confusion operations along vertical direction
    for(int j = BYTE_MATRIX_COLUMNS - 1; j >= 0; j--){
        uint16_t shiftDistance = shiftDistanceSequence[-- shiftDistanceSequenceIndex] % byteMatrixRows;

        for(int i = 0; i < BYTE_MATRIX_ROWS; i++)
            byteMatrix[((i - shiftDistance) + byteMatrixRows) % byteMatrixRows][j] = tempByteMatrix[i][j];
    }
    cloneByteMatrix(byteMatrix, tempByteMatrix);
    
    //perform inverse confusion operations along horizontal direction
    for(int i = BYTE_MATRIX_ROWS - 1; i >= 0; i--){
        uint16_t shiftDistance = shiftDistanceSequence[-- shiftDistanceSequenceIndex] % byteMatrixColumns;

        for(int j = 0; j < BYTE_MATRIX_COLUMNS; j++)
            byteMatrix[i][((j - shiftDistance) + byteMatrixColumns) % byteMatrixColumns] = tempByteMatrix[i][j];
    }
    cloneByteMatrix(byteMatrix, tempByteMatrix);

    return shiftDistanceSequenceIndex;
}

//display the byte matrix
void displayByteMatrix(uint8_t byteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS]){
    for(int i = 0; i < BYTE_MATRIX_ROWS; i++){
        for(int j = 0; j < BYTE_MATRIX_COLUMNS; j++)
            printf("%-3d ", byteMatrix[i][j]);

        printf("\n");
    }
}

//convert iteration results into byte sequence
void convertIterationResultsToByteSequence(double * iterationResultArray1, double * iterationResultArray2, uint8_t * byteSequence, int numberOfIterationResults){
    uint8_t uint8tDataArray1[numberOfIterationResults * NUMBER_OF_UINT8T_DATA_EXTRACTED_FROM_EACH_ITERATION_RESULT];
    uint8_t uint8tDataArray2[numberOfIterationResults * NUMBER_OF_UINT8T_DATA_EXTRACTED_FROM_EACH_ITERATION_RESULT];

    //convert iteration results into uint8_t data
    for(int i = 0; i < numberOfIterationResults; i++){
        uint8_t *p = &uint8tDataArray1[i * NUMBER_OF_UINT8T_DATA_EXTRACTED_FROM_EACH_ITERATION_RESULT];
        memcpy(p, &iterationResultArray1[i], NUMBER_OF_UINT8T_DATA_EXTRACTED_FROM_EACH_ITERATION_RESULT);

        p = &uint8tDataArray2[i * NUMBER_OF_UINT8T_DATA_EXTRACTED_FROM_EACH_ITERATION_RESULT];
        memcpy(p, &iterationResultArray2[i], NUMBER_OF_UINT8T_DATA_EXTRACTED_FROM_EACH_ITERATION_RESULT);
    }

    //generate byte sequence
    for(int i = 0; i < numberOfIterationResults * NUMBER_OF_UINT8T_DATA_EXTRACTED_FROM_EACH_ITERATION_RESULT; i++)
        byteSequence[i] = uint8tDataArray1[i] ^ uint8tDataArray2[i];
}

//generate diffusion seed
uint8_t generateDiffusionSeed(double * initialCondition1, double controlParameter1, double * initialCondition2, double controlParameter2){
    //generate iteration result
    double iterationResult1 = PLCM(* initialCondition1, controlParameter1);
    double iterationResult2 = PLCM(* initialCondition2, controlParameter2);

    uint8_t uint8tData1, uint8tData2;

    memcpy(&uint8tData1, &iterationResult1, 1);
    memcpy(&uint8tData2, &iterationResult2, 1);

    return uint8tData1 ^ uint8tData2;
}

//perform diffusion operations on the byte matrix
int diffusion(uint8_t byteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS], uint8_t tempByteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS], uint8_t * byteSequence, uint8_t diffusionSeed, int byteSequenceIndex){
    int prei, prej;

    for(int i = 0; i < BYTE_MATRIX_ROWS; i++)
        for(int j = 0; j < BYTE_MATRIX_COLUMNS; j++){
            if(i == 0 && j == 0){
                byteMatrix[i][j] = byteSequence[byteSequenceIndex] ^ ((tempByteMatrix[i][j] + byteSequence[byteSequenceIndex]) % 256) ^ diffusionSeed;
                byteSequenceIndex ++;
            }
            
            else if(i != 0 && j == 0){
                prei = i - 1;
                prej = BYTE_MATRIX_COLUMNS - 1;

                byteMatrix[i][j] = byteSequence[byteSequenceIndex] ^ ((tempByteMatrix[i][j] + byteSequence[byteSequenceIndex]) % 256) ^ byteMatrix[prei][prej];
                byteSequenceIndex ++;
            }   

            else{
                prei = i;
                prej = j - 1;

                byteMatrix[i][j] = byteSequence[byteSequenceIndex] ^ ((tempByteMatrix[i][j] + byteSequence[byteSequenceIndex]) % 256) ^ byteMatrix[prei][prej];
                byteSequenceIndex ++;
            }
        }
    
    cloneByteMatrix(byteMatrix, tempByteMatrix);

    return byteSequenceIndex;
}

//perform inverse diffusion operations on the byte matrix
int inverseDiffusion(uint8_t byteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS], uint8_t tempByteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS], uint8_t * byteSequence, uint8_t diffusionSeed, int byteSequenceIndex){
    int prei, prej;

    for(int i = BYTE_MATRIX_ROWS - 1; i >= 0; i--)
        for(int j = BYTE_MATRIX_COLUMNS - 1; j >= 0; j--){
            if(i == 0 && j == 0){
                byteSequenceIndex = byteSequenceIndex - 1;
                byteMatrix[i][j] = ((byteSequence[byteSequenceIndex] ^ tempByteMatrix[i][j] ^ diffusionSeed) + 256 - byteSequence[byteSequenceIndex]);
            }
            
            else if(i != 0 && j == 0){
                prei = i - 1;
                prej = BYTE_MATRIX_COLUMNS - 1;

                byteSequenceIndex = byteSequenceIndex - 1;
                byteMatrix[i][j] = ((byteSequence[byteSequenceIndex] ^ tempByteMatrix[i][j] ^ tempByteMatrix[prei][prej]) + 256 - byteSequence[byteSequenceIndex]);
            }

            else{
                prei = i;
                prej = j - 1;

                byteSequenceIndex = byteSequenceIndex - 1;
                byteMatrix[i][j] = ((byteSequence[byteSequenceIndex] ^ tempByteMatrix[i][j] ^ tempByteMatrix[prei][prej]) + 256 - byteSequence[byteSequenceIndex]);
            }
        }

    cloneByteMatrix(byteMatrix, tempByteMatrix);

    return byteSequenceIndex;
}

//dispaly SHA-256 Hash
void dispalySHA256Hash(uint8_t * SHA256HashResultsArray){
    int index = 0;

    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 8; j++)
            printf("0x%02X ", SHA256HashResultsArray[index ++]);

        printf("\n");
    }
}

//draw waveform of the audio samples 
void drawWaveform(SDL_Renderer * renderer, uint8_t * sampleBuffer, int bufferSize, int position) {
    switch(position) {
        case 0: SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255); break;   
        case 1: SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255); break;   
        case 2: SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255); break;    
        default: SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);     
    }

    int sectionHeight = WINDOW_HEIGHT / 3;
    int midY = (position * sectionHeight) + (sectionHeight / 2);  
    int amplitude = sectionHeight / 2 - 5; 

    for (int i = 0; i < WINDOW_WIDTH; i++) {
        int sampleIndex = (i * bufferSize) / WINDOW_WIDTH;
        int16_t sample = *((int16_t*)(sampleBuffer + sampleIndex));
        int y = midY + (sample * amplitude) / 32768;

        SDL_Rect point = {i - 1, y - 1, 3, 3}; 
        SDL_RenderFillRect(renderer, &point);
    }
}

void drawByteMatrix(Mat frame, Mat canvas, int position) {

    int sectionHeight = canvas.rows / 3;
    int y_start = position * sectionHeight;

    double scale = (double)sectionHeight / frame.rows;
    Mat resizedImg;
    resize(frame, resizedImg, Size(), scale, scale);

    int x_start = (canvas.cols - resizedImg.cols) / 2;

    Mat roi = canvas(Rect(x_start, y_start, resizedImg.cols, resizedImg.rows));
    resizedImg.copyTo(roi);

    if (position > 0) {
        line(canvas, 
             Point(0, position * sectionHeight), 
             Point(canvas.cols, position * sectionHeight), 
             Scalar(0, 0, 0), 2);
    }
}
