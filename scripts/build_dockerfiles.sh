#!/usr/bin/env bash
set -euo pipefail

push=false
if [[ "${1:-}" == --push ]]; then
  push=true
  shift
fi
if [[ $# -lt 2 ]]; then
  echo "Usage: $0 [--push] REGISTRY RELEASE_TAG [IMAGE ...]" >&2
  exit 2
fi
registry="${1%/}"
tag="$2"
shift 2
if [[ -z "$registry" || ! "$tag" =~ ^[[:alnum:]_][[:alnum:]_.-]{0,127}$ ]]; then
  echo "Registry must be nonempty and RELEASE_TAG must be a valid Docker tag" >&2
  exit 2
fi
repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo"
if [[ $# -eq 0 ]]; then
  for dockerfile in dockerfiles/*.Dockerfile; do
    name="${dockerfile##*/}"
    set -- "$@" "${name%.Dockerfile}"
  done
fi
for name in "$@"; do
  if [[ ! "$name" =~ ^[a-z0-9_]+$ || ! -f "dockerfiles/$name.Dockerfile" ]]; then
    echo "Unknown image: $name" >&2
    exit 2
  fi
done
for name in "$@"; do
  image="$registry/$name:$tag"
  docker build --platform linux/amd64 -t "$image" -f "dockerfiles/$name.Dockerfile" .
  if "$push"; then
    docker push "$image"
    docker image inspect --format '{{join .RepoDigests "\n"}}' "$image"
  fi
done
