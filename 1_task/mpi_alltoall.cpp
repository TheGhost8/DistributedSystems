#include <iostream>
#include <vector>
#include <mpi.h>

/*
0  - 1  - 2  - 3
|    |    |    |
4  - 5  - 6  - 7
|    |    |    |
8  - 9  - 10 - 11
|    |    |    |
12 - 13 - 14 - 15
*/

static size_t MatrixRank = 16;


void addRecievedData(std::vector<int> &data, std::vector<int> &buffer) {
    if (data.size() != MatrixRank) {
        std::cout << "Data size != MatrixRank" << std::endl;
        return;
    }
    if (data.size() != buffer.size()) {
        std::cout << "Data size != buffer size" << std::endl;
        return;
    }

    for (size_t i = 0; i < data.size(); ++i) {
        if (buffer[i] != 0) {
            data[i] = buffer[i];
        }
    }
}

int dataSum(std::vector<int> &data)
{
    int sum = 0;
    
    for (auto &val : data) {
        sum += val;
    }

    return sum;
}

void sendMPIData(size_t thread_to_send, std::vector<int> &data) {
    MPI_Request request;
    MPI_Isend(data.data(), data.size(), MPI_INT, thread_to_send, 0, MPI_COMM_WORLD, &request);
}

void recieveMPIData(size_t thread_to_recv, std::vector<int> &data, size_t proc_num) {
    MPI_Status status;
    MPI_Request request;
    std::vector<int> buffer;
    buffer.resize(MatrixRank);

    MPI_Recv(&buffer[0], MatrixRank, MPI_INT, thread_to_recv, 0, MPI_COMM_WORLD, &status);
    // MPI_Irecv(&buffer[0], MatrixRank, MPI_INT, thread_to_recv, 0, MPI_COMM_WORLD, &request);
    if (status.MPI_ERROR != MPI_SUCCESS) {
        std::cout << "Something went wrong during MPI_Recv from " + std::to_string(thread_to_recv) + " on " + std::to_string(proc_num) + " proc" << std::endl;
        return;
    }
    addRecievedData(data, buffer);
}

void createMPIData(std::vector<std::vector<int>> &all_data, size_t thread_num, std::vector<size_t> &proc_threads)
{
    all_data.resize(MatrixRank);

    for (size_t i = 0; i < all_data.size(); ++i) {
        all_data[i].resize(MatrixRank);
        for (size_t j = 0; j < all_data[i].size(); ++j) {
            if ((proc_threads[i] == thread_num) && (j == i)) {
                all_data[i][j] = i + 1;
            }
            else {
                all_data[i][j] = 0;
            }
        }
    }
}

void checkPrint(std::vector<std::vector<int>> &all_data)
{
    for (auto &proc_data : all_data) {
        std::cout << std::to_string(dataSum(proc_data)) + " ";
    }
    std::cout << std::endl;
}

void firstStep(size_t thread_num, std::vector<size_t> &proc_threads, std::vector<std::vector<int>> &all_data) {
    size_t proc_num_send = 0, proc_num_recv = 0, p = 0;

    for (p = 0; p < 4; ++p) {
        if (thread_num == proc_threads[p]) {
            proc_num_send = p + 4;
            sendMPIData(proc_threads[proc_num_send], all_data[p]);
        }
    }
    for (p = 12; p < 16; ++p) {
        if (thread_num == proc_threads[p]) {
            proc_num_send = p - 4;
            sendMPIData(proc_threads[proc_num_send], all_data[p]);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    for (p = 4; p < 12; ++p) {
        if (thread_num == proc_threads[p]) {
            if (p < 8) {
                proc_num_recv = p - 4;    
            }
            else {
                proc_num_recv = p + 4;
            }
            recieveMPIData(proc_threads[proc_num_recv], all_data[p], p);
        }
    }
}

void secondStep(size_t thread_num, std::vector<size_t> &proc_threads, std::vector<std::vector<int>> &all_data) {
    size_t proc_num_send = 0, proc_num_recv = 0, p = 0;

    for (p = 4; p <= 8; p += 4) {
        if (thread_num == proc_threads[p]) {
            proc_num_send = p + 1;
            sendMPIData(proc_threads[proc_num_send], all_data[p]);
        }
        if (thread_num == proc_threads[p + 3]) {
            proc_num_send = p + 2;
            sendMPIData(proc_threads[proc_num_send], all_data[p + 3]);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    for (p = 5; p <= 9; p += 4) {
        if (thread_num == proc_threads[p]) {
            proc_num_recv = p - 1;
            recieveMPIData(proc_threads[proc_num_recv], all_data[p], p);
        }
        if (thread_num == proc_threads[p + 1]) {
            proc_num_recv = p + 2;
            recieveMPIData(proc_threads[proc_num_recv], all_data[p + 1], p + 1);
        }
    }
}

void thirdStep(size_t thread_num, std::vector<size_t> &proc_threads, std::vector<std::vector<int>> &all_data) {
    size_t proc_num_send = 0, proc_num_recv = 0, p = 0;

    if (thread_num == proc_threads[5]) {
        proc_num_send = 9;
        sendMPIData(proc_threads[proc_num_send], all_data[5]);
    }
    if (thread_num == proc_threads[10]) {
        proc_num_send = 6;
        sendMPIData(proc_threads[proc_num_send], all_data[10]);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (thread_num == proc_threads[9]) {
        proc_num_recv = 5;
        recieveMPIData(proc_threads[proc_num_recv], all_data[9], 9);
    }
    if (thread_num == proc_threads[6]) {
        proc_num_recv = 10;
        recieveMPIData(proc_threads[proc_num_recv], all_data[6], 6);
    }

    if (thread_num == proc_threads[9]) {
        proc_num_send = 5;
        sendMPIData(proc_threads[proc_num_send], all_data[9]);
    }
    if (thread_num == proc_threads[6]) {
        proc_num_send = 10;
        sendMPIData(proc_threads[proc_num_send], all_data[6]);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (thread_num == proc_threads[5]) {
        proc_num_recv = 9;
        recieveMPIData(proc_threads[proc_num_recv], all_data[5], 5);
    }
    if (thread_num == proc_threads[10]) {
        proc_num_recv = 6;
        recieveMPIData(proc_threads[proc_num_recv], all_data[10], 10);
    }
}

void fourthStep(size_t thread_num, std::vector<size_t> &proc_threads, std::vector<std::vector<int>> &all_data) {
    size_t proc_num_send = 0, proc_num_recv = 0, p = 0;

    if (thread_num == proc_threads[6]) {
        proc_num_send = 5;
        sendMPIData(proc_threads[proc_num_send], all_data[6]);
    }
    if (thread_num == proc_threads[9]) {
        proc_num_send = 10;
        sendMPIData(proc_threads[proc_num_send], all_data[9]);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (thread_num == proc_threads[5]) {
        proc_num_recv = 6;
        recieveMPIData(proc_threads[proc_num_recv], all_data[5], 5);
    }
    if (thread_num == proc_threads[10]) {
        proc_num_recv = 9;
        recieveMPIData(proc_threads[proc_num_recv], all_data[10], 10);
    }

    if (thread_num == proc_threads[5]) {
        proc_num_send = 6;
        sendMPIData(proc_threads[proc_num_send], all_data[5]);
    }
    if (thread_num == proc_threads[10]) {
        proc_num_send = 9;
        sendMPIData(proc_threads[proc_num_send], all_data[10]);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (thread_num == proc_threads[6]) {
        proc_num_recv = 5;
        recieveMPIData(proc_threads[proc_num_recv], all_data[6], 6);
    }
    if (thread_num == proc_threads[9]) {
        proc_num_recv = 10;
        recieveMPIData(proc_threads[proc_num_recv], all_data[9], 9);
    }
}

void fifthStep(size_t thread_num, std::vector<size_t> &proc_threads, std::vector<std::vector<int>> &all_data) {
    size_t proc_num_send = 0, proc_num_recv = 0, p = 0;

    if (thread_num == proc_threads[5]) {
        proc_num_send = 4;
        sendMPIData(proc_threads[proc_num_send], all_data[5]);
        proc_num_send = 1;
        sendMPIData(proc_threads[proc_num_send], all_data[5]);
    }
    if (thread_num == proc_threads[6]) {
        proc_num_send = 2;
        sendMPIData(proc_threads[proc_num_send], all_data[6]);
        proc_num_send = 7;
        sendMPIData(proc_threads[proc_num_send], all_data[6]);
    }
    if (thread_num == proc_threads[9]) {
        proc_num_send = 8;
        sendMPIData(proc_threads[proc_num_send], all_data[9]);
        proc_num_send = 13;
        sendMPIData(proc_threads[proc_num_send], all_data[9]);
    }
    if (thread_num == proc_threads[10]) {
        proc_num_send = 11;
        sendMPIData(proc_threads[proc_num_send], all_data[10]);
        proc_num_send = 14;
        sendMPIData(proc_threads[proc_num_send], all_data[10]);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (thread_num == proc_threads[1]) {
        proc_num_recv = 5;
        recieveMPIData(proc_threads[proc_num_recv], all_data[1], 1);
    }
    if (thread_num == proc_threads[4]) {
        proc_num_recv = 5;
        recieveMPIData(proc_threads[proc_num_recv], all_data[4], 4);
    }
    if (thread_num == proc_threads[2]) {
        proc_num_recv = 6;
        recieveMPIData(proc_threads[proc_num_recv], all_data[2], 2);
    }
    if (thread_num == proc_threads[7]) {
        proc_num_recv = 6;
        recieveMPIData(proc_threads[proc_num_recv], all_data[7], 7);
    }
    if (thread_num == proc_threads[8]) {
        proc_num_recv = 9;
        recieveMPIData(proc_threads[proc_num_recv], all_data[8], 8);
    }
    if (thread_num == proc_threads[13]) {
        proc_num_recv = 9;
        recieveMPIData(proc_threads[proc_num_recv], all_data[13], 13);
    }
    if (thread_num == proc_threads[11]) {
        proc_num_recv = 10;
        recieveMPIData(proc_threads[proc_num_recv], all_data[11], 11);
    }
    if (thread_num == proc_threads[14]) {
        proc_num_recv = 10;
        recieveMPIData(proc_threads[proc_num_recv], all_data[14], 14);
    }
}

void sixthStep(size_t thread_num, std::vector<size_t> &proc_threads, std::vector<std::vector<int>> &all_data) {
    size_t proc_num_send = 0, proc_num_recv = 0, p = 0;

    for (p = 1; p <= 13; p += 12) {
        if (thread_num == proc_threads[p]) {
            proc_num_send = p - 1;
            sendMPIData(proc_threads[proc_num_send], all_data[p]);
        }
        if (thread_num == proc_threads[p + 1]) {
            proc_num_send = p + 2;
            sendMPIData(proc_threads[proc_num_send], all_data[p + 1]);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    for (p = 0; p <= 12; p += 12) {
        if (thread_num == proc_threads[p]) {
            proc_num_recv = p + 1;
            recieveMPIData(proc_threads[proc_num_recv], all_data[p], p);
        }
        if (thread_num == proc_threads[p + 3]) {
            proc_num_recv = p + 2;
            recieveMPIData(proc_threads[proc_num_recv], all_data[p + 3], p + 3);
        }
    }
}

void testMatrix(size_t thread_num, std::vector<size_t> &proc_threads) {
    std::vector<std::vector<int>> all_data;

    createMPIData(all_data, thread_num, proc_threads);

    firstStep(thread_num, proc_threads, all_data);

    secondStep(thread_num, proc_threads, all_data);

    thirdStep(thread_num, proc_threads, all_data);

    fourthStep(thread_num, proc_threads, all_data);

    fifthStep(thread_num, proc_threads, all_data);

    sixthStep(thread_num, proc_threads, all_data);

    if (thread_num == 0) {
        checkPrint(all_data);
    }
}

int main(int argc, char** argv)
{
    MPI_Init(NULL, NULL);
    double timeStart = 0, timeFinish = 0;
    int rank = 0, size = 0;

    MPI_Comm_size (MPI_COMM_WORLD, &size);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    std::vector<size_t> thread_for_proc;

    for (size_t i = 0; i < MatrixRank; ++i) {
        proc_threads.push_back(std::make_pair(i, i % size));
        thread_for_proc.push_back(i % size);
    }

    timeStart = MPI_Wtime();
    testMatrix(rank, thread_for_proc);
    timeFinish = MPI_Wtime();

    if (rank == 0) {
        std::cout << "Time elapsed: " << std::to_string(timeFinish - timeStart) << std::endl;
    }

    MPI_Finalize();
    return 0;
}
