
!
!----
!

  subroutine send_c(sendbuf, sendcount, dest, sendtag) !! FS FS ASKI

  use my_mpi

  implicit none

  integer :: dest,sendtag,sendcount,ier
  complex, dimension(sendcount):: sendbuf

  call MPI_SEND(sendbuf,sendcount,MPI_COMPLEX,dest,sendtag, &
       my_local_mpi_comm_world,ier)

  end subroutine send_c

!
!----
!

  subroutine recv_c(recvbuf, recvcount, source, recvtag) !! FS FS ASKI

  use my_mpi

  implicit none

  integer :: source,recvtag,recvcount,ier
  complex, dimension(recvcount):: recvbuf

  call MPI_RECV(recvbuf,recvcount,MPI_COMPLEX,source,recvtag, &
       my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_c

