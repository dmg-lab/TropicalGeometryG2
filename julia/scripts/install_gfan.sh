#!/bin/bash
# Build gfan 0.8 beta without root, from the source and bug-fix patches of the
# Fedora 44 package. The compiler, GMP, cddlib and TBB come from conda-forge
# through micromamba. Installs into ~/opt/gfan and links gfan and
# gfan_tropicalprevariety into ~/bin; PREFIX, BINDIR and JOBS override these.
#   bash scripts/install_gfan.sh
set -euo pipefail

PREFIX=${PREFIX:-$HOME/opt/gfan}
BINDIR=${BINDIR:-$HOME/bin}
JOBS=${JOBS:-16}
DEPS=$PREFIX/deps
SOURCE=https://math.au.dk/~jensen/software/gfan/gfan0.8beta.tar.gz
SOURCE_SHA256=fa7884e5f317c50f8fb4f37bcf5d419f0fd5f7b90d6037349d1957ea73cebbee
FEDORA=https://src.fedoraproject.org/rpms/gfan/raw/8316b24c39f962ea9a438f0db79f595d37fff796/f
# The shared-library, SoPlex and C++20 patches only concern Fedora's packaging.
PATCHES=(
    "warning 48f79ed18378a4ceb9eb0314257644358f32fd48821e23633fb9bde3cd84bd19"
    "gcd 05d505b8e54c236ad5a6b5139919cbd39c1ead45b39c8ce50b3175dd1b980874"
    "multiplicities c802a096cf4dd46e01fb155f487c510603dda079c8ea0c2449e446740ada45c7"
    "empty-coneStack fcfe61ba0f017777738a3dea535e69fc07985e9b28ae6ad95904c98950a2f607"
)

fetch() {
    curl -fsSL "$1" -o "$2"
    echo "$3  $2" | sha256sum -c --quiet
}

mkdir -p "$PREFIX/build" "$PREFIX/bin" "$BINDIR"
cd "$PREFIX/build"

if [ ! -x "$PREFIX/micromamba" ]; then
    curl -fsSL https://micro.mamba.pm/api/micromamba/linux-64/latest |
        tar -xj -C "$PREFIX" --strip-components=1 bin/micromamba
fi
if [ ! -d "$DEPS" ]; then
    MAMBA_ROOT_PREFIX=$PREFIX/mamba "$PREFIX/micromamba" create -y -q -p "$DEPS" \
        -c conda-forge --override-channels gcc=14 gxx=14 make patch gmp cddlib tbb-devel
fi
export PATH=$DEPS/bin:$PATH

rm -rf gfan0.8beta
fetch "$SOURCE" gfan0.8beta.tar.gz "$SOURCE_SHA256"
tar -xzf gfan0.8beta.tar.gz
for entry in "${PATCHES[@]}"; do
    read -r name sha256 <<< "$entry"
    fetch "$FEDORA/gfan-$name.patch" "gfan-$name.patch" "$sha256"
    patch -s -p1 -d gfan0.8beta < "gfan-$name.patch"
done

# Fedora's optimization flags; -march=native, so build on the machine that runs it.
make -s -C gfan0.8beta -j"$JOBS" CC=gcc CXX=g++ cddnoprefix=true \
    OPTFLAGS="-O2 -march=native -DGMPRATIONAL -DNDEBUG" \
    CDD_INCLUDEOPTIONS="-I$DEPS/include/cddlib" CDD_LINKOPTIONS="-L$DEPS/lib -lcddgmp" \
    GMP_INCLUDEOPTIONS="-I$DEPS/include" GMP_LINKOPTIONS="-L$DEPS/lib -lgmp" \
    TBB_INCLUDEOPTIONS="-I$DEPS/include" \
    PLATFORM_LINKOPTIONS="-rdynamic -L$DEPS/lib -ltbb -Wl,-rpath,$DEPS/lib"

install -m 755 gfan0.8beta/gfan "$PREFIX/bin/gfan"
ln -sf "$PREFIX/bin/gfan" "$BINDIR/gfan"
ln -sf "$PREFIX/bin/gfan" "$BINDIR/gfan_tropicalprevariety"
usage=$("$BINDIR/gfan_tropicalprevariety" --help 2>&1 || true)
[[ $usage == *--halfopenrestrictions* ]] || { echo "the installed gfan is not 0.8" >&2; exit 1; }
echo "gfan 0.8 beta installed: $BINDIR/gfan_tropicalprevariety"
