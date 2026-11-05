find ~ -type d -name ".git" -prune -print0 2>/dev/null |
while IFS= read -r -d '' gitdir; do
    repo="${gitdir%/.git}"

    git -C "$repo" stash list --format="%gd" |
    while read -r stash; do

        # Search tracked files in the stash
        result=$(git -C "$repo" grep -n -I -w "numB" "$stash" 2>/dev/null)

        if [ -n "$result" ]; then
            echo ""
            echo "======================================"
            echo "REPOSITORY: $repo"
            echo "STASH: $stash (tracked files)"
            echo "$result"
        fi

        # Search untracked files in the stash (third parent)
        if git -C "$repo" rev-parse "$stash^3" >/dev/null 2>&1; then
            result=$(git -C "$repo" grep -n -I -w "numB" "$stash^3" 2>/dev/null)

            if [ -n "$result" ]; then
                echo ""
                echo "======================================"
                echo "REPOSITORY: $repo"
                echo "STASH: $stash (untracked files)"
                echo "$result"
            fi
        fi

    done
done