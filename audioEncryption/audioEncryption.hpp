#ifndef __AUDIOENCRYPTION_HPP__
#define __AUDIOENCRYPTION_HPP__

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <SDL2/SDL.h>
#include <pulse/error.h>
#include <pulse/simple.h>
#include <openssl/sha.h>
#include "opencv2/opencv.hpp"

using namespace cv;

//macros for initializing capture and playback devices
#define CHANNELS                                                    2
#define SAMPLE_WIDTH                                                2
#define SAMPLING_RATE                                               44100

//macros for encrypting and decrypting the samples
#define FRAMES_IN_EACH_BLOCK                                        1024
#define SAMPLE_BUFFER_SIZE                                          FRAMES_IN_EACH_BLOCK * CHANNELS * SAMPLE_WIDTH
#define BYTE_MATRIX_ROWS                                            32
#define BYTE_MATRIX_COLUMNS                                         SAMPLE_BUFFER_SIZE / BYTE_MATRIX_ROWS

#define PRE_ITERATIONS                                              200
#define ROUNDS_OF_CONFUSION_OPERATIONS                              3
#define ROUNDS_OF_DIFFUSION_CONFUSION_OPERATIONS                    3
#define NUMBER_OF_UINT16T_DATA_EXTRACTED_FROM_EACH_ITERATION_RESULT 3
#define NUMBER_OF_UINT8T_DATA_EXTRACTED_FROM_EACH_ITERATION_RESULT  6

//macros for SDL
#define WINDOW_WIDTH                                                800
#define WINDOW_HEIGHT                                               600

//get cpu time
double getCPUSecond(void);

//convert the linear samples into byte matrix
void convertSampleBufferToByteMatrix(uint8_t * sampleBuffer, uint8_t byteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS]);

//convert the byte matrix into the linear samples
void convertByteMatrixToSampleBuffer(uint8_t byteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS], uint8_t * sampleBuffer);

//calculate the SHA-256 hash of the captured samples
void calculateSHA256Hash(uint8_t * sampleBuffer, uint8_t * SHA256HashResultsArray);

//reconsturct initial conditions and control parameters using SHA-256 hash of the captured samples
double reconstructParameter(double parameter, uint8_t * SHA256HashResultsArray);

//iterate the PLCM and return the iteration result
double PLCM(double initialCondition, double controlParameter);

//convert iteration results into shift distance sequence
void convertIterationResultsToShiftDistanceSequence(double * iterationResultArray1, double * iterationResultArray2, uint16_t * shiftDistanceSequence, int numberOfIterationResults);

//clone byte matrix
void cloneByteMatrix(uint8_t sourceByteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS], uint8_t targetBytematrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS]);

//perform confusion operations on the byte matrix
int confusion(uint8_t byteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS], uint8_t tempByteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS], uint16_t * shiftDistanceSequence, int shiftDistanceSequenceIndex);

//perform inverse confusion operations on the byte matrix
int inverseConfusion(uint8_t byteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS], uint8_t tempByteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS], uint16_t * shiftDistanceSequence, int shiftDistanceSequenceIndex);

//display the byte matrix
void displayByteMatrix(uint8_t byteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS]);

//convert iteration results into byte sequence
void convertIterationResultsToByteSequence(double * iterationResultArray1, double * iterationResultArray2, uint8_t * byteSequence, int numberOfIterationResults);

//generate diffusion seed
uint8_t generateDiffusionSeed(double * initialCondition1, double controlParameter1, double * initialCondition2, double controlParameter2);

//perform diffusion operations on the byte matrix
int diffusion(uint8_t byteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS], uint8_t tempByteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS], uint8_t * byteSequence, uint8_t diffusionSeed, int byteSequenceIndex);

//perform inverse diffusion operations on the byte matrix
int inverseDiffusion(uint8_t byteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS], uint8_t tempByteMatrix[BYTE_MATRIX_ROWS][BYTE_MATRIX_COLUMNS], uint8_t * byteSequence, uint8_t diffusionSeed, int byteSequenceIndex);

//dispaly SHA-256 Hash
void dispalySHA256Hash(uint8_t * SHA256HashResultsArray);

//draw waveform of the audio samples 
void drawWaveform(SDL_Renderer * renderer, uint8_t * sampleBuffer, int bufferSize, int position);

//draw the figure of byte matrix
void drawByteMatrix(Mat frame, Mat canvas, int position);

#endif
