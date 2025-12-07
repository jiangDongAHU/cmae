g++ audioEncryption.cpp -o out -lpulse-simple -lpulse -lssl -lcrypto -lSDL2 -Wno-deprecated-declarations `pkg-config --cflags --libs opencv4`
