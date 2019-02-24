function max_samples = find_length_of_notes(notes)
max_samples = 0;
    for ii = 1:length(notes)
        if notes{ii}.start + notes{ii}.duration > max_samples
            max_samples = notes{ii}.start + notes{ii}.duration;
        end
    end
end