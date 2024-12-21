# chmod +x deploy.sh
echo 'removing data, objects and executables'
rm *.dat *.o *.x

echo 'generating objects and executables'
make module_algebra.o
make module_subroutines.o
make von_karman_swirling_flow.o
make von_karman_swirling_flow.x

echo 'running'
./von_karman_swirling_flow.x