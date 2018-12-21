#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "core.h"
#include "mpi.h"

#define EXTERNAL -1
#define MAX_NEIGHBOURS 7
#define NO_TAG 0

#define FAIL_EXCHANGE 3

extern int nproc;
extern int myrank;
extern int i_am_the_master;
extern MPI_Comm comm_grid;

static double *sendbuf;
int *message_len;
int **messages;

void initialize(CSRMatrix A, int *Rows, int *Part, int part_len) {

  ni = A.N;
  
  int *Glob2Loc = (int *)malloc(part_len * sizeof(int));
  for (int i = 0; i < part_len; i++) {
    Glob2Loc[i] = EXTERNAL;
  }

  // fill local rows
  for (int i = 0; i < ni; i++) {
    Glob2Loc[Rows[i]] = i;
  }

  // if (i_am_the_master) {
  //   printf("\n\n-----------Proc %d-------------\n\n", myrank);
  //   printf("Glob2Loc after filling local:\n[");
  //   for (int i = 0; i < part_len; i++) {
  //     printf("%d, ", Glob2Loc[i]);
  //   }
  //   printf("]\n");
  // }

  // length of messages
  message_len = (int *)malloc(nproc * sizeof(int));
  for (int i = 0; i < nproc; i++) {
    message_len[i] = 0;
  }

  // calculating messages length
  nh = 0;
  for (int i = 0; i < A.N; i++) {
    for (int j = A.IA[i]; j < A.IA[i+1]; j++) {
      if (Glob2Loc[A.JA[j]] == EXTERNAL) {
        int owner = Part[A.JA[j]];
        message_len[owner]++;
        nh++;
      }
    }
  }
  // if (i_am_the_master) {
  //   printf("ni=%d, nh=%d\n", ni, nh);
  // }
  sendbuf = (double *)malloc(nh * sizeof(double));

  // number of elements currently in each message
  int *msg_size = (int *)malloc(nproc * sizeof(int));
  for (int i = 0; i < nproc; i++) {
    msg_size[i] = 0;
  }

  // lists of elements to be exchanged with each process
  messages = (int **)malloc(nproc * sizeof(int *));
  for (int i = 0; i < nproc; i++) {
    if (message_len[i] == 0) {
      messages[i] = NULL;
    } else {
      messages[i] = (int *)malloc(message_len[i] * sizeof(int));
    }
  }

  // calculating elements to receive
  for (int i = 0; i < A.N; i++) {
    for (int j = A.IA[i]; j < A.IA[i+1]; j++) {
      int glob_elem = A.JA[j];
      if (Glob2Loc[glob_elem] == EXTERNAL) {
        int owner = Part[glob_elem];
        int k = msg_size[owner];
        messages[owner][k] = glob_elem;
        msg_size[owner]++;
      }
    }
  }

  // int token;
  // if (myrank == 0) {
  //   printf("Receive messages [proc %d]:\n", myrank);
  //   for (int i = 0; i < nproc; i++) {
  //     if (message_len[i] > 0) {
  //       printf("Proc %d [%d]: [", i, message_len[i]);
  //       for (int j = 0; j < message_len[i]; j++) {
  //         printf("%d, ", messages[i][j]);
  //       }
  //       printf("]\n");
  //     }
  //   }
  //   MPI_Send(&token, 1, MPI_INT, myrank+1, 0, comm_grid);
  // } else {
  //   MPI_Recv(&token, 1, MPI_INT, myrank-1, 0, comm_grid, MPI_STATUS_IGNORE);
  //   printf("Receive messages [proc %d]:\n", myrank);
  //   for (int i = 0; i < nproc; i++) {
  //     if (message_len[i] > 0) {
  //       printf("Proc %d [%d]: [", i, message_len[i]);
  //       for (int j = 0; j < message_len[i]; j++) {
  //         printf("%d, ", messages[i][j]);
  //       }
  //       printf("]\n");
  //     }
  //   }
  //   if (myrank != nproc-1) {
  //     MPI_Send(&token, 1, MPI_INT, myrank+1, 0, comm_grid);
  //   }
  // }
  // MPI_Barrier(comm_grid);

  // assert that messages are completed
  for (int i = 0; i < nproc; i++) {
    assert(msg_size[i] == message_len[i]);
    assert(msg_size[i] > 0 || messages[i] == NULL);
  }

  // map elements to receive into Glob2Loc
  int idx = ni; // halo indices start at ni
  for (int i = 0; i < nproc; i++) {
    for (int j = 0; j < message_len[i]; j++) {
      Glob2Loc[messages[i][j]] = idx;
      idx++;
    }
  }

  // calculating messages to send
  for (int i = 0; i < nproc; i++) {
    msg_size[i] = 0;
  }
  for (int i = 0; i < A.N; i++) {
    for (int j = A.IA[i]; j < A.IA[i+1]; j++) {
      int glob_elem = A.JA[j];
      // if neighbour is in halo
      if (Glob2Loc[glob_elem] >= A.N) {
        int owner = Part[glob_elem];
        int k = msg_size[owner];
        // then we send current row
        messages[owner][k] = i;
        msg_size[owner]++;
      }
    }
  }

  // if (myrank == 0) {
  //   printf("\n\n");
  //   printf("Send messages [proc %d]:\n", myrank);
  //   for (int i = 0; i < nproc; i++) {
  //     if (message_len[i] > 0) {
  //       printf("Proc %d [%d]: [", i, message_len[i]);
  //       for (int j = 0; j < message_len[i]; j++) {
  //         printf("%d, ", Rows[messages[i][j]]);
  //       }
  //       printf("]\n");
  //     }
  //   }
  //   MPI_Send(&token, 1, MPI_INT, myrank+1, 0, comm_grid);
  // } else {
  //   MPI_Recv(&token, 1, MPI_INT, myrank-1, 0, comm_grid, MPI_STATUS_IGNORE);
  //   printf("Send messages [proc %d]:\n", myrank);
  //   for (int i = 0; i < nproc; i++) {
  //     if (message_len[i] > 0) {
  //       printf("Proc %d [%d]: [", i, message_len[i]);
  //       for (int j = 0; j < message_len[i]; j++) {
  //         printf("%d, ", Rows[messages[i][j]]);
  //       }
  //       printf("]\n");
  //     }
  //   }
  //   if (myrank != nproc-1) {
  //     MPI_Send(&token, 1, MPI_INT, myrank+1, 0, comm_grid);
  //   }
  // }
  // MPI_Barrier(comm_grid);

  // asserting that messages are complete
  for (int i = 0; i < nproc; i++) {
    assert(msg_size[i] == message_len[i]);
    assert(msg_size[i] > 0 || messages[i] == NULL);
  }

  // mapping JA to local numeration
  for (int i = 0; i < A.N; i++) {
    for (int j = A.IA[i]; j < A.IA[i+1]; j++) {
      int glob_elem = A.JA[j];
      int loc_elem = Glob2Loc[glob_elem];
      A.JA[j] = loc_elem;
    }
  }

  free(msg_size);
  free(Glob2Loc);
}

void do_exchange(Vector V) {

  MPI_Request rq[MAX_NEIGHBOURS * 2];
  MPI_Status st[MAX_NEIGHBOURS * 2];

  int start = 0;
  int req_no = 0;
  for (int i = 0; i < nproc; i++) {
    if (message_len[i] == 0) {
      continue;
    }
    // writing elements into sendbuf
    for (int j = 0; j < message_len[i]; j++) {
      int loc_elem = messages[i][j];
      sendbuf[start + j] = V.data[loc_elem];
    }
    MPI_Isend(sendbuf + start, message_len[i], MPI_DOUBLE, i, NO_TAG, comm_grid, rq + req_no);
    start += message_len[i];
    req_no++;
  }
  assert(V.nh == start);

  start = V.ni;  // halo starts from ni
  for (int i = 0; i < nproc; i++) {
    if (message_len[i] == 0) {
      continue;
    }
    // receive elements into the vector
    MPI_Irecv(V.data + start, message_len[i], MPI_DOUBLE, i, NO_TAG, comm_grid, rq + req_no);
    start += message_len[i];
    req_no++;
  }

  // waiting for all communications to finish
  int retval = MPI_Waitall(req_no, (MPI_Request *)rq, (MPI_Status *)st);
  if (retval == MPI_ERR_IN_STATUS) {
    fprintf(stderr, "%s\n", "One of the exchanges was completed with an error");
    exit(FAIL_EXCHANGE);
  }
}

void finalize() {
  free(sendbuf);
  free(message_len);
  for (int i = 0; i < nproc; i++) {
    if (messages[i] != NULL) {
      free(messages[i]);
    }
  }
}