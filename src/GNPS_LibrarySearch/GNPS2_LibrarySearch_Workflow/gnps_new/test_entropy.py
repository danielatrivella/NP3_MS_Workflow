import pytest
import numpy as np
from .entropy import entropy_similarity

def test_entropy_similarity_standard():
    peaks1 = np.array([[50, 8.0], [60, 100.0], [80, 50.0], [100, 50.0]], dtype=np.float32)
    peaks2 = np.array([[55, 38.0], [70, 50.0], [80, 66.0], [90, 100.0]], dtype=np.float32)
    score, n_matches = entropy_similarity(peaks1, peaks2, tolerance=0.05, sqrt_transform=True, penalty=0)
    print(f"Standard Score: {score:.3f}, Matches: {n_matches}")
    assert score == pytest.approx(0.26056383550167084)
    assert n_matches == 1

def test_entropy_similarity_penalty():
    peaks1 = np.array([[50, 8.0], [60, 100.0], [80, 50.0], [100, 50.0]], dtype=np.float32)
    peaks2 = np.array([[55, 38.0], [70, 50.0], [80, 66.0], [90, 100.0]], dtype=np.float32)
    score, n_matches = entropy_similarity(peaks1, peaks2, tolerance=0.05, sqrt_transform=True, penalty=0.6)
    print(f"Penalty Score: {score:.3f}, Matches: {n_matches}")
    assert score == pytest.approx(0.3425292372703552)
    assert n_matches == 1

def test_entropy_similarity_with_shift():
    peaks1 = np.array([[50, 8.0], [60, 100.0], [80, 50.0], [100, 50.0]], dtype=np.float32)
    peaks2 = np.array([[55, 38.0], [70, 50.0], [80, 66.0], [90, 100.0]], dtype=np.float32)
    score, n_matches = entropy_similarity(peaks1, peaks2, tolerance=0.05, sqrt_transform=True, penalty=0.6, shift=10)
    print(f"Enhanced Reverse Score with Shift: {score:.3f}, Matches: {n_matches}")
    assert score == pytest.approx(0.6542002707719803)
    assert n_matches == 2

def test_entropy_similarity_empty_spectra():
    peaks1 = np.array([], dtype=np.float32).reshape(0, 2)
    peaks2 = np.array([], dtype=np.float32).reshape(0, 2)
    score, n_matches = entropy_similarity(peaks1, peaks2, tolerance=0.05, sqrt_transform=True, penalty=0)
    assert score == 0.0
    assert n_matches == 0

def test_entropy_similarity_unsorted_spectra():
    unsorted_peaks1 = np.array([[100, 50.0], [80, 50.0], [60, 100.0], [50, 8.0]], dtype=np.float32)
    unsorted_peaks2 = np.array([[90, 100.0], [80, 66.0], [70, 50.0], [55, 38.0]], dtype=np.float32)
    unsorted_score, unsorted_n_matches = entropy_similarity(unsorted_peaks1, unsorted_peaks2, tolerance=0.05, sqrt_transform=True, penalty=0, require_sorted=False)
    sorted_peaks1   = np.array([[50, 8.0], [60, 100.0], [80, 50.0], [100, 50.0]], dtype=np.float32)
    sorted_peaks2   = np.array([[55, 38.0], [70, 50.0], [80, 66.0], [90, 100.0]], dtype=np.float32)
    sorted_score, sorted_n_matches = entropy_similarity(sorted_peaks1, sorted_peaks2, tolerance=0.05, sqrt_transform=True, penalty=0, require_sorted=True)
    print(f"Unsorted Score: {unsorted_score:.3f}, Matches: {unsorted_n_matches}")
    print(f"Sorted Score: {sorted_score:.3f}, Matches: {sorted_n_matches}")
    assert pytest.approx(unsorted_score, 0.2499243052576287) == sorted_score
    assert unsorted_n_matches == sorted_n_matches
    assert unsorted_score == sorted_score

def test_entropy_similarity_sorted_requirement():
    peaks1 = np.array([[100, 50.0], [80, 50.0], [60, 100.0], [50, 8.0]], dtype=np.float32)
    peaks2 = np.array([[90, 100.0], [80, 66.0], [70, 50.0], [55, 38.0]], dtype=np.float32)
    with pytest.raises(ValueError):
        entropy_similarity(peaks1, peaks2, tolerance=0.05, sqrt_transform=True, penalty=0, require_sorted=True)
