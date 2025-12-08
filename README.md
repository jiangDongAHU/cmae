## Multi-User Secure Audio Communication via Video Encryption Framework: A Cross-Modal Cryptographic Approach

### Jingyi Ni, Zean Feng, Jian Zhang, Zijian Cui, Hongkai Tian, and Dong Jiang

#### Develoment enviroment

* CPU             : Intel Xeon Gold 6226R @ 2.9GHz
* Memory          : 32GB
* Operating system: Ubuntu 22.04
* Pipewire        : 15.99.0
* OpenCV          : 4.5.4

#### File description

* audioEncryption: This directory contains the source code implementing the encryption and decryption procedures of the proposed protocol. Executing the setup.sh script within this directory automatically compiles the source code and initiates the program.
* demo.webm      : This demonstration video showcases the program’s execution workflow. The program acquires audio samples from the default microphone, encrypts them using the proposed protocol, subsequently decrypts the encrypted data, and plays back the recovered audio. In parallel, it displays the amplitudes of the original, encrypted, and decrypted audio samples, and visualizes the corresponding byte matrices.
