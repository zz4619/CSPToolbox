! Read-only extraction of CP2's compiled space-group table for this benchmark.
program dump_cp2_groups
use crystal_symmetry
implicit none
integer i,j,k
do i=1,NSpaceSupported
  write(*,'(a,3i8,1x,a)') 'GROUP ',i,SpaceSupported(i)%Number,SpaceSupported(i)%NSymOp,trim(SpaceSupported(i)%Name)
  write(*,'(3i8,3es25.16)') SpaceSupported(i)%len1,SpaceSupported(i)%len2,SpaceSupported(i)%len3,&
       SpaceSupported(i)%ang1,SpaceSupported(i)%ang2,SpaceSupported(i)%ang3
  do j=1,SpaceSupported(i)%NSymOp
    write(*,'(12es25.16)') SpaceSupported(i)%SymRot(:,:,j),SpaceSupported(i)%SymTrn(:,j)
  enddo
enddo
end program
